/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <math.h>
#include "core/encseq.h"
#include "core/codetype.h"
#include "core/arraydef.h"
#include "core/showtime.h"
#include "core/mathsupport.h"
#include "core/radix-intsort.h"
#include "sfx-suffixer.h"
#include "firstcodes.h"

typedef struct
{
  unsigned long code, *ptr;
} GtIndexwithcode;

GT_DECLAREARRAYSTRUCT(GtIndexwithcode);

typedef struct
{
  unsigned int units, shiftright;
  unsigned long entries;
  GtCodetype mask;
  GtBitsequence *bits;
} Gtmarksubstring;

static void gt_marksubstring_init(Gtmarksubstring *mark,unsigned int numofchars,
                                  unsigned int kmersize, bool usesuffix,
                                  unsigned int units)
{
  gt_assert(kmersize >= units);
  mark->units = units;
  mark->entries = gt_power_for_small_exponents(numofchars,units);
  gt_assert(kmersize >= units + units);
  if (usesuffix)
  {
    mark->shiftright = 0;
    mark->mask = (GtCodetype) (mark->entries-1);
  } else
  {
    mark->shiftright = GT_MULT2(kmersize - units);
    mark->mask = ~(GtCodetype) 0;
  }
  GT_INITBITTAB(mark->bits,mark->entries);
}

static void gt_marksubstring_mark(Gtmarksubstring *mark,GtCodetype code)
{
  code = (code >> (GtCodetype) mark->shiftright) & mark->mask;

  gt_assert(code < mark->entries);
  if (!GT_ISIBITSET(mark->bits,code))
  {
    GT_SETIBIT(mark->bits,code);
  }
}

static bool gt_marksubstring_checkmark(const Gtmarksubstring *mark,
                                       GtCodetype code)
{
  code = (code >> (GtCodetype) mark->shiftright) & mark->mask;

  return GT_ISIBITSET(mark->bits,code) ? true : false;
}

typedef struct
{
  unsigned long differentcodes,
                firstcodehits,
                countsequences,
                numofsequences,
                binsearchcache_unknownwidth,
                binsearchcache_unknowncount,
                *codebuffer,
                codebuffer_allocated,
                codebuffer_nextfree,
                codebuffer_total,
                *allfirstcodes,
                *countocc;
  GtUlongPair *codeposbuffer;
  GtArrayGtIndexwithcode binsearchcache;
  unsigned int binsearchcache_depth,
               flushcount;
  Gtmarksubstring markprefix,
                  marksuffix;
  unsigned long *tempforradixsort;
  uint32_t *suftabseqnum;
  uint8_t *suftaboffset;
} GtFirstcodesinfo;

static void gt_storefirstcodes(void *processinfo,
                               GT_UNUSED bool firstinrange,
                               GT_UNUSED unsigned long pos,
                               GtCodetype code)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) processinfo;

  gt_assert(firstinrange);
  gt_assert(fci->allfirstcodes != NULL &&
            fci->countsequences < fci->numofsequences);
  fci->allfirstcodes[fci->countsequences] = code;
  fci->countsequences++;
}

#ifdef SKDEBUG
static void checkcodesorder(const unsigned long *tab,unsigned long len,
                            bool allowequal)
{
  unsigned long idx;

  for (idx=1UL; idx < len; idx++)
  {
    gt_assert(tab[idx-1] < tab[idx] || (allowequal && tab[idx-1] == tab[idx]));
  }
}
#endif

static unsigned long gt_remdups_in_sorted_array(GtFirstcodesinfo *fci)
{
  unsigned long *storeptr, *readptr;

  if (fci->numofsequences > 0)
  {
    unsigned long numofdifferentcodes;

    fci->countocc = gt_calloc((size_t) fci->numofsequences,
                              sizeof (*fci->countocc));
    fci->countocc[0] = 1UL;
    gt_marksubstring_mark(&fci->markprefix,fci->allfirstcodes[0]);
    gt_marksubstring_mark(&fci->marksuffix,fci->allfirstcodes[0]);
    for (storeptr = fci->allfirstcodes, readptr = fci->allfirstcodes+1;
         readptr < fci->allfirstcodes + fci->numofsequences;
         readptr++)
    {
      if (*storeptr != *readptr)
      {
        storeptr++;

        *storeptr = *readptr;
      }
      fci->countocc[(unsigned long) (storeptr - fci->allfirstcodes)]++;
      gt_marksubstring_mark(&fci->markprefix,*readptr);
      gt_marksubstring_mark(&fci->marksuffix,*readptr);
    }
    numofdifferentcodes = (unsigned long) (storeptr - fci->allfirstcodes + 1);
    if (numofdifferentcodes < fci->numofsequences)
    {
      /* reduce the memory requirement, as the duplicated elements are not
         needed */
      fci->allfirstcodes = gt_realloc(fci->allfirstcodes,
                                      sizeof (*fci->allfirstcodes) *
                                      numofdifferentcodes);
      fci->countocc = gt_realloc(fci->countocc,
                                 sizeof (*fci->countocc) *
                                 (numofdifferentcodes+1));
#ifdef SKDEBUG
      checkcodesorder(fci->allfirstcodes,numofdifferentcodes,false);
#endif
    }
    gt_assert(fci->countocc != NULL);
    return numofdifferentcodes;
  }
  return 0;
}

static void gt_firstcodes_halves_rek(GtFirstcodesinfo *fci,
                                     unsigned long left,unsigned long right,
                                     unsigned int depth,unsigned int maxdepth)
{
  unsigned long mid;

  gt_assert(left <= right);
  mid = left + GT_DIV2(right-left);
  if (depth < maxdepth)
  {
    gt_firstcodes_halves_rek(fci,left,mid-1,depth+1,maxdepth);
  } else
  {
    fci->binsearchcache_unknownwidth += (unsigned long) (mid-left);
    fci->binsearchcache_unknowncount++;
  }
  gt_assert(fci->binsearchcache.nextfreeGtIndexwithcode <
            fci->binsearchcache.allocatedGtIndexwithcode);
  fci->binsearchcache.spaceGtIndexwithcode[fci->binsearchcache.
                                           nextfreeGtIndexwithcode]
                                           .ptr = fci->allfirstcodes + mid;
  fci->binsearchcache.spaceGtIndexwithcode[fci->binsearchcache.
                                           nextfreeGtIndexwithcode++]
                                           .code = fci->allfirstcodes[mid];
  if (depth < maxdepth)
  {
    gt_firstcodes_halves_rek(fci,mid+1,right,depth+1,maxdepth);
  } else
  {
    fci->binsearchcache_unknownwidth += (unsigned long) (right-mid);
    fci->binsearchcache_unknowncount++;
  }
}

static void gt_firstcodes_halves(GtFirstcodesinfo *fci,
                                 unsigned int maxdepth)
{
  fci->binsearchcache.nextfreeGtIndexwithcode = 0;
  fci->binsearchcache.allocatedGtIndexwithcode = (1UL << (maxdepth+1)) - 1;
  if (fci->binsearchcache.allocatedGtIndexwithcode < fci->differentcodes)
  {
    size_t allocbytes
      = sizeof (*fci->binsearchcache.spaceGtIndexwithcode)
                * fci->binsearchcache.allocatedGtIndexwithcode;
    printf("size of binsearch cache: %lu\n",(unsigned long) allocbytes);
    fci->binsearchcache.spaceGtIndexwithcode = gt_malloc(allocbytes);
    gt_assert(fci->differentcodes > 0);
    gt_firstcodes_halves_rek(fci,0,fci->differentcodes - 1,0,maxdepth);
    gt_assert(fci->binsearchcache.nextfreeGtIndexwithcode
              == fci->binsearchcache.allocatedGtIndexwithcode);
    printf("average size of uncached range: %.2f\n",
            (double) fci->binsearchcache_unknownwidth/
            fci->binsearchcache_unknowncount);
#undef FIRSTCODEDEBUG
#ifdef FIRSTCODEDEBUG
    {
      unsigned long idx;

      for (idx=0; idx < fci->binsearchcache.nextfreeGtIndexwithcode;
           idx++)
      {
        printf("%lu %lu\n",
             (unsigned long)
             (fci->binsearchcache.spaceGtIndexwithcode[idx].ptr -
             fci->allfirstcodes),
             fci->binsearchcache.spaceGtIndexwithcode[idx].code);
      }
    }
#endif
  }
}

const unsigned long *gt_firstcodes_find(const GtFirstcodesinfo *fci,
                                        unsigned long code)
{
  const unsigned long *leftptr = NULL, *midptr, *rightptr = NULL;
  unsigned int depth;

  if (fci->binsearchcache.spaceGtIndexwithcode != NULL)
  {
    const GtIndexwithcode *leftic, *midic, *rightic;

    leftic = fci->binsearchcache.spaceGtIndexwithcode;
    rightic = fci->binsearchcache.spaceGtIndexwithcode +
              fci->binsearchcache.nextfreeGtIndexwithcode - 1;
    for (depth = 0; /* Nothing */; depth++)
    {
      midic = leftic + GT_DIV2((unsigned long) (rightic-leftic));
      if (code < midic->code)
      {
        if (depth < fci->binsearchcache_depth)
        {
          rightic = midic - 1;
        } else
        {
          gt_assert(leftic->ptr != NULL && rightic->ptr != NULL);
          if (leftic > fci->binsearchcache.spaceGtIndexwithcode)
          {
            leftptr = (leftic-1)->ptr + 1;
          } else
          {
            leftptr = fci->allfirstcodes;
          }
          rightptr = rightic->ptr - 1;
          break;
        }
      } else
      {
        if (code > midic->code)
        {
          if (depth < fci->binsearchcache_depth)
          {
            leftic = midic + 1;
          } else
          {
            gt_assert(leftic->ptr != NULL && rightic->ptr != NULL);
            leftptr = leftic->ptr + 1;
            if (rightic < fci->binsearchcache.spaceGtIndexwithcode +
                     fci->binsearchcache.nextfreeGtIndexwithcode-1)
            {
              rightptr = (rightic+1)->ptr - 1;
            } else
            {
              rightptr = fci->allfirstcodes +
                         fci->differentcodes - 1;
            }
            break;
          }
        } else
        {
          return midic->ptr;
        }
      }
    }
    gt_assert(leftptr != NULL && rightptr != NULL);
  } else
  {
    depth = 0;
    leftptr = fci->allfirstcodes;
    rightptr = fci->allfirstcodes + fci->differentcodes - 1;
  }
  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    if (code < *midptr)
    {
      rightptr = midptr - 1;
    } else
    {
      if (code > *midptr)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr;
      }
    }
    depth++;
  }
  return NULL;
}

static unsigned long gt_firstcodes_accumulatecounts_merge(
                                        GtFirstcodesinfo *fci,
                                        const unsigned long *querystream_fst,
                                        const unsigned long *subjectstream_fst)
{
  unsigned long found = 0;
  const unsigned long *query = querystream_fst,
                      *subject = subjectstream_fst,
                      *querystream_lst = fci->codebuffer
                                         + fci->codebuffer_nextfree,
                      *subjectstream_lst = fci->allfirstcodes
                                           + fci->differentcodes;

  while (query < querystream_lst && subject < subjectstream_lst)
  {
    if (*query < *subject)
    {
      query++;
    } else
    {
      if (*query > *subject)
      {
        subject++;
      } else
      {
        fci->countocc[(unsigned long)
                                 (subject - fci->allfirstcodes)]++;
        query++;
        found++;
      }
    }
  }
  return found;
}

static void gt_firstcodes_accumulatecounts_flush(GtFirstcodesinfo
                                                   *fci)
{
  const unsigned long *ptr;
  unsigned long *vptr;

  gt_assert(fci->allfirstcodes != NULL);
  gt_radixsort_GtUlong_linear(true,fci->codebuffer,
                              fci->tempforradixsort,fci->codebuffer_nextfree);
#ifdef SKDEBUG
  checkcodesorder(fci->codebuffer,fci->codebuffer_nextfree,
                  true);
#endif
  fci->codebuffer_total += fci->codebuffer_nextfree;
  for (vptr = fci->codebuffer;
       vptr < fci->codebuffer +
              fci->codebuffer_nextfree;
       vptr++)
  {
    ptr = gt_firstcodes_find(fci,*vptr);
    if (ptr != NULL)
    {
      fci->firstcodehits += gt_firstcodes_accumulatecounts_merge(fci,vptr, ptr);
      break;
    }
  }
  printf("%u ",fci->flushcount++);
  (void) fflush(stdout);
  fci->codebuffer_nextfree = 0;
}

static void gt_firstcodes_accumulatecounts(void *processinfo,
                                           GT_UNUSED bool firstinrange,
                                           GT_UNUSED unsigned long pos,
                                           GtCodetype code)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) processinfo;

  if (fci->codebuffer_nextfree  == fci->codebuffer_allocated)
  {
    gt_firstcodes_accumulatecounts_flush(fci);
  }
  if (gt_marksubstring_checkmark(&fci->markprefix,code) &&
      gt_marksubstring_checkmark(&fci->marksuffix,code))
  {
    gt_assert (fci->codebuffer_nextfree <
               fci->codebuffer_allocated);
    fci->codebuffer[fci->codebuffer_nextfree++] = code;
  }
}

#ifdef INSERTION
static unsigned long gt_firstcodes_insertsuffix_merge(
                                        GtFirstcodesinfo *fci,
                                        const unsigned long *querystream_fst,
                                        const unsigned long *subjectstream_fst)
{
  unsigned long found = 0, idx;
  const unsigned long *query = querystream_fst,
                      *subject = subjectstream_fst,
                      *querystream_lst = fci->codebuffer
                          + fci->codebuffer_nextfree,
                      *subjectstream_lst
                        = fci->allfirstcodes +
                          fci->differentcodes;

  while (query < querystream_lst && subject < subjectstream_lst)
  {
    if (*query < *subject)
    {
      query++;
    } else
    {
      if (*query > *subject)
      {
        subject++;
      } else
      {
        idx = --fci->countocc[(unsigned long)
                               (subject - fci->allfirstcodes)];
        query++;
        found++;
      }
    }
  }
  return found;
}
#endif

static void storefirstcodes_partialsum(GtFirstcodesinfo *fci)
{
  unsigned long idx;

  for (idx = 1UL; idx < fci->differentcodes; idx++)
  {
    fci->countocc[idx] += fci->countocc[idx-1];
  }
  fci->countocc[fci->differentcodes] = fci->countocc[fci->differentcodes-1];
}

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME)                        firstcodes_##NAME
#define firstcodes_ARRAY_GET(ARR,RELIDX)       ARR[RELIDX]
#define firstcodes_ARRAY_SET(ARR,RELIDX,VALUE) ARR[RELIDX] = VALUE

typedef unsigned long QSORTNAME(Sorttype);

#include "match/qsort-direct.gen"

void storefirstcodes_getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                                   unsigned int kmersize)
{
  GtTimer *timer = NULL;
  GtFirstcodesinfo fci;
  size_t sizeforcodestable, numofsuftabentries;
  unsigned int numofchars = gt_encseq_alphabetnumofchars(encseq);
  const unsigned int markprefixunits = 14U;

  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description("insert first codes into "
                                                   "array");
    gt_timer_start(timer);
  }
  fci.numofsequences = gt_encseq_num_of_sequences(encseq);
  gt_assert(numofchars == 4U);
  sizeforcodestable = sizeof (*fci.allfirstcodes) *
                      fci.numofsequences;
  printf("# use array of size %lu\n",(unsigned long) sizeforcodestable);
  fci.allfirstcodes = gt_malloc(sizeforcodestable);
  fci.differentcodes = 0;
  fci.countsequences = 0;
  fci.firstcodehits = 0;
  fci.binsearchcache_unknownwidth = 0;
  fci.binsearchcache_unknowncount = 0;
  GT_INITARRAY(&fci.binsearchcache,GtIndexwithcode);
  getencseqkmers_twobitencoding(encseq,
                                GT_READMODE_FORWARD,
                                kmersize,
                                kmersize,
                                true,
                                gt_storefirstcodes,
                                &fci,
                                NULL,
                                NULL);
  gt_assert(fci.numofsequences == fci.countsequences);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "sorting the codes",stdout);
  }
  QSORTNAME(gt_direct_qsort)
             (6UL, false,
             fci.allfirstcodes,fci.numofsequences);
  gt_marksubstring_init(&fci.markprefix,numofchars,kmersize,false,
                        markprefixunits);
  gt_marksubstring_init(&fci.marksuffix,numofchars,kmersize,
                        true,markprefixunits);
  fci.differentcodes = gt_remdups_in_sorted_array(&fci);
  printf("# number of different codes=%lu (%.4f) in %lu sequences\n",
          fci.differentcodes,
          (double) fci.differentcodes/fci.numofsequences,
          fci.countsequences);
  fci.binsearchcache_depth = 15U;
  fci.flushcount = 0;
  fci.codebuffer_total = 0;
  gt_firstcodes_halves(&fci,fci.binsearchcache_depth);
  fci.codebuffer_allocated = 3000000UL;
  fci.codebuffer_nextfree = 0;
  fci.codebuffer = gt_malloc(sizeof (*fci.codebuffer)
                             * fci.codebuffer_allocated);
  gt_assert(fci.codebuffer_allocated > 0);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "accumulate counts",stdout);
  }
  fci.tempforradixsort = gt_malloc((size_t) fci.codebuffer_allocated
                                   * sizeof (*fci.tempforradixsort));
  getencseqkmers_twobitencoding(encseq,
                                GT_READMODE_FORWARD,
                                kmersize,
                                45U,
                                false,
                                gt_firstcodes_accumulatecounts,
                                &fci,
                                NULL,
                                NULL);
  gt_firstcodes_accumulatecounts_flush(&fci);
  printf("\ncodebuffer_total=%lu (%.2f)\n",
          fci.codebuffer_total,
          (double) fci.codebuffer_total/
                   gt_encseq_total_length(encseq));
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "compute partial sums",stdout);
  }
  printf("# firstcodehits=%lu (%.2f)\n",fci.firstcodehits,
                                      (double) fci.firstcodehits/
                                      gt_encseq_total_length(encseq));
  storefirstcodes_partialsum(&fci);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "insert suffixes",stdout);
  }
  numofsuftabentries = (size_t) fci.firstcodehits;
  gt_free(fci.tempforradixsort);
  gt_free(fci.codebuffer);
  fci.firstcodehits = 0;
  gt_assert(fci.codebuffer_nextfree == 0);
  fci.codebuffer_total = 0;
  fci.codebuffer_allocated /= 2UL;
  fci.codeposbuffer = gt_malloc(sizeof (*fci.codeposbuffer)
                                * fci.codebuffer_allocated);
  fci.suftabseqnum = gt_malloc(numofsuftabentries * sizeof (*fci.suftabseqnum));
  fci.suftaboffset = gt_malloc(numofsuftabentries * sizeof (*fci.suftaboffset));
  /*
  getencseqkmers_twobitencoding(encseq,
                                GT_READMODE_FORWARD,
                                kmersize,
                                45U,
                                false,
                                gt_firstcodes_insertsuffixes,
                                &fci,
                                NULL,
                                NULL);
  */
  GT_FREEARRAY(&fci.binsearchcache,GtIndexwithcode);
  gt_free(fci.codeposbuffer);
  gt_free(fci.allfirstcodes);
  gt_free(fci.countocc);
  gt_free(fci.markprefix.bits);
  gt_free(fci.marksuffix.bits);
  gt_free(fci.suftabseqnum);
  gt_free(fci.suftaboffset);
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
}
