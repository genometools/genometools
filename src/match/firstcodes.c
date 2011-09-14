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
#include "core/log_api.h"
#include "core/spacepeak.h"
#include "core/spacecalc.h"
#include "core/fa.h"
#include "sfx-suffixer.h"
#include "sfx-shortreadsort.h"
#include "spmsuftab.h"
#include "esa-spmsk.h"
#include "firstcodes.h"
#include "sfx-maprange.h"
#include "stamp.h"

typedef struct
{
  unsigned long code, *ptr;
} GtIndexwithcode;

GT_DECLAREARRAYSTRUCT(GtIndexwithcode);

typedef struct
{
  size_t size;
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
  mark->size = sizeof (GtBitsequence) * GT_NUMOFINTSFORBITS(mark->entries);
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
  unsigned long firstcodehits,
                firstcodeposhits,
                countsequences,
                numofsequences,
                binsearchcache_unknownwidth,
                binsearchcache_unknowncount,
                maxbucketsize,
                *codebuffer,
                codebuffer_allocated,
                codebuffer_nextfree,
                codebuffer_total;
  GtUlongPair *codeposbuffer,
              *tempcodeposforradixsort;
  GtArrayGtIndexwithcode binsearchcache;
  unsigned int binsearchcache_depth,
               flushcount;
  Gtmarksubstring markprefix,
                  marksuffix;
  unsigned long *tempcodeforradixsort;
  GtSpmsuftab *spmsuftab;
  size_t workspace;
  size_t size_to_split;
  unsigned long differentcodes,
                *allfirstcodes,
                *countocc;
  GtSfxmappedrange *mappedcountocc,
                   *mappedallfirstcodes,
                   *mappedmarkprefix;
} GtFirstcodesinfo;

/* call the following function after computing the partial sums */

unsigned long gt_firstcodes_get_leftborder(const GtFirstcodesinfo *fci,
                                           unsigned long idx)
{
  gt_assert(idx <= fci->differentcodes);
  return fci->countocc[idx];
}

size_t gt_firstcodes_size_to_split(const GtFirstcodesinfo *fci)
{
  return fci->size_to_split;
}

unsigned long gt_firstcodes_numofallcodes(const GtFirstcodesinfo *fci)
{
  return fci->differentcodes;
}

unsigned long gt_firstcodes_findfirstlarger(const GtFirstcodesinfo *fci,
                                            unsigned long suftaboffset)
{
  unsigned long left = 0, right = fci->differentcodes, mid, midval,
                found = fci->differentcodes;

  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_firstcodes_get_leftborder(fci,mid);
    if (suftaboffset == midval)
    {
      return mid;
    }
    if (suftaboffset < midval)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  return found;
}

unsigned long gt_firstcodes_mapped_range_size(const GtFirstcodesinfo *fci,
                                              unsigned long minindex,
                                              unsigned long maxindex)
{
  size_t idx;
  GtSfxmappedrange *maptab[3];
  unsigned long sumsize = 0;

  maptab[0] = fci->mappedcountocc;
  maptab[1] = fci->mappedallfirstcodes;
  maptab[2] = fci->mappedmarkprefix;;
  for (idx = 0; idx < sizeof (maptab)/sizeof (maptab[0]); idx++)
  {
    if (maptab[idx] == NULL)
    {
      return (unsigned long) fci->size_to_split;
    }
    sumsize += gt_Sfxmappedrange_size_mapped(fci->mappedcountocc,minindex,
                                             maxindex);
  }
  return sumsize;
}

unsigned long gt_firstcodes_idx_code(const GtFirstcodesinfo *fci,
                                     unsigned long idx)
{
  gt_assert(idx < fci->differentcodes);
  return fci->allfirstcodes[idx];
}

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

    fci->countocc = gt_calloc((size_t) (fci->numofsequences+1),
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
      fci->size_to_split
        += sizeof (*fci->allfirstcodes) * numofdifferentcodes
           + sizeof (*fci->countocc) * (numofdifferentcodes+1);
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

static size_t gt_firstcodes_halves(GtFirstcodesinfo *fci,
                                   unsigned int maxdepth)
{
  fci->binsearchcache.nextfreeGtIndexwithcode = 0;
  fci->binsearchcache.allocatedGtIndexwithcode = (1UL << (maxdepth+1)) - 1;
  if (fci->binsearchcache.allocatedGtIndexwithcode < fci->differentcodes)
  {
    size_t allocbytes
      = sizeof (*fci->binsearchcache.spaceGtIndexwithcode)
                * fci->binsearchcache.allocatedGtIndexwithcode;
    fci->binsearchcache.spaceGtIndexwithcode = gt_malloc(allocbytes);
    gt_assert(fci->differentcodes > 0);
    gt_firstcodes_halves_rek(fci,0,fci->differentcodes - 1,0,maxdepth);
    gt_assert(fci->binsearchcache.nextfreeGtIndexwithcode
              == fci->binsearchcache.allocatedGtIndexwithcode);
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
    return allocbytes;
  }
  return 0;
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
        fci->countocc[(unsigned long) (subject - fci->allfirstcodes)]++;
        query++;
        found++;
      }
    }
  }
  return found;
}

static void gt_firstcodes_accumulatecounts_flush(GtFirstcodesinfo *fci)
{
  const unsigned long *ptr;
  unsigned long *vptr;

  gt_assert(fci->allfirstcodes != NULL);
  gt_radixsort_GtUlong_linear(false,fci->codebuffer,
                              fci->tempcodeforradixsort,
                              fci->codebuffer_nextfree);
#ifdef SKDEBUG
  checkcodesorder(fci->codebuffer,fci->codebuffer_nextfree,
                  true);
#endif
  fci->codebuffer_total += fci->codebuffer_nextfree;
  for (vptr = fci->codebuffer;
       vptr < fci->codebuffer + fci->codebuffer_nextfree;
       vptr++)
  {
    ptr = gt_firstcodes_find(fci,*vptr);
    if (ptr != NULL)
    {
      fci->firstcodehits += gt_firstcodes_accumulatecounts_merge(fci,vptr, ptr);
      break;
    }
  }
#ifdef FLUSHCOUNT
  if (fci->flushcount == 0)
  {
    printf("# ");
  }
  printf("%u ",fci->flushcount++);
  (void) fflush(stdout);
#endif
  fci->codebuffer_nextfree = 0;
}

static void gt_firstcodes_accumulatecounts(void *processinfo,
                                           bool firstinrange,
                                           GT_UNUSED unsigned long pos,
                                           GtCodetype code)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) processinfo;

  if (!firstinrange &&
      gt_marksubstring_checkmark(&fci->markprefix,code) &&
      gt_marksubstring_checkmark(&fci->marksuffix,code))
  {
    if (fci->codebuffer_nextfree == fci->codebuffer_allocated)
    {
      gt_firstcodes_accumulatecounts_flush(fci);
    }
    gt_assert (fci->codebuffer_nextfree < fci->codebuffer_allocated);
    fci->codebuffer[fci->codebuffer_nextfree++] = code;
  }
}

static unsigned long gt_firstcodes_insertsuffixes_merge(
                                        GtFirstcodesinfo *fci,
                                        const GtUlongPair *querystream_fst,
                                        const unsigned long *subjectstream_fst)
{
  unsigned long found = 0, idx;
  const GtUlongPair *query = querystream_fst,
                      *querystream_lst = fci->codeposbuffer +
                                         fci->codebuffer_nextfree;
  const unsigned long *subject = subjectstream_fst,
                      *subjectstream_lst = fci->allfirstcodes +
                                           fci->differentcodes;

  while (query < querystream_lst && subject < subjectstream_lst)
  {
    if (query->a < *subject)
    {
      query++;
    } else
    {
      if (query->a > *subject)
      {
        subject++;
      } else
      {
        idx = --fci->countocc[(unsigned long) (subject - fci->allfirstcodes)];
        gt_assert(idx < fci->firstcodehits + fci->numofsequences);
        gt_spmsuftab_set(fci->spmsuftab,idx,query->b);
        query++;
        found++;
      }
    }
  }
  return found;
}

static void gt_firstcodes_insertsuffixes_flush(GtFirstcodesinfo *fci)
{
  const unsigned long *ptr;
  GtUlongPair *vptr;

  gt_assert(fci->allfirstcodes != NULL);
  gt_radixsort_GtUlongPair_linear(false,fci->codeposbuffer,
                                  fci->tempcodeposforradixsort,
                                  fci->codebuffer_nextfree);
#ifdef SKDEBUG
  checkcodeposorder(fci->codeposbuffer,fci->codebuffer_nextfree,true);
#endif
  fci->codebuffer_total += fci->codebuffer_nextfree;
  for (vptr = fci->codeposbuffer;
       vptr < fci->codeposbuffer +
              fci->codebuffer_nextfree;
       vptr++)
  {
    ptr = gt_firstcodes_find(fci,vptr->a);
    if (ptr != NULL)
    {
      fci->firstcodeposhits += gt_firstcodes_insertsuffixes_merge(fci,vptr,ptr);
      break;
    }
  }
#ifdef FLUSHCOUNT
  if (fci->flushcount == 0)
  {
    printf("# ");
  }
  printf("%u ",fci->flushcount++);
  (void) fflush(stdout);
#endif
  fci->codebuffer_nextfree = 0;
}

static void gt_firstcodes_insertsuffixes(void *processinfo,
                                         GT_UNUSED bool firstinrange,
                                         unsigned long pos,
                                         GtCodetype code)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) processinfo;

  if (gt_marksubstring_checkmark(&fci->markprefix,code) &&
      gt_marksubstring_checkmark(&fci->marksuffix,code))
  {
    if (fci->codebuffer_nextfree == fci->codebuffer_allocated)
    {
      gt_firstcodes_insertsuffixes_flush(fci);
    }
    gt_assert (fci->codebuffer_nextfree < fci->codebuffer_allocated);
    fci->codeposbuffer[fci->codebuffer_nextfree].a = code;
    fci->codeposbuffer[fci->codebuffer_nextfree++].b = pos;
  }
}

static void storefirstcodes_partialsum(GtFirstcodesinfo *fci)
{
  unsigned long idx;

  fci->maxbucketsize = fci->countocc[0];
  for (idx = 1UL; idx < fci->differentcodes; idx++)
  {
    if (fci->maxbucketsize < fci->countocc[idx])
    {
      fci->maxbucketsize = fci->countocc[idx];
    }
    fci->countocc[idx] += fci->countocc[idx-1];
  }
  fci->countocc[fci->differentcodes] = fci->countocc[fci->differentcodes-1];
}

static void gt_firstcodes_checksuftab_bucket(const GtEncseq *encseq,
                                             GtReadmode readmode,
                                             GtEncseqReader *esr1,
                                             GtEncseqReader *esr2,
                                             unsigned long previous,
                                             bool previousdefined,
                                             const unsigned long *suftabbuffer,
                                             GT_UNUSED const uint16_t
                                               *lcptab_bucket,
                                             unsigned long numberofsuffixes)
{
  unsigned long idx, current, maxlcp,
                totallength = gt_encseq_total_length(encseq);
  int cmp;
  const bool specialsareequal = false, specialsareequalatdepth0 = false;
  const unsigned long depth = 0;

  gt_assert(!previousdefined || previous < totallength);
  for (idx = 0; idx < numberofsuffixes; idx++)
  {
    current = suftabbuffer[idx];
    if (previousdefined && idx < totallength)
    {
      gt_assert(current < totallength);
      cmp = gt_encseq_check_comparetwosuffixes(encseq,
                                               readmode,
                                               &maxlcp,
                                               specialsareequal,
                                               specialsareequalatdepth0,
                                               depth,
                                               previous,
                                               current,
                                               esr1,
                                               esr2);
      gt_assert(cmp <= 0);
      gt_assert(idx == 0 || maxlcp == (unsigned long) lcptab_bucket[idx]);
    }
    previous = current;
    previousdefined = true;
  }
}

static void gt_firstcodes_sortremaining(const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        GtSpmsuftab *spmsuftab,
                                        unsigned long maxbucketsize,
                                        const unsigned long *countocc,
                                        unsigned long differentcodes,
                                        unsigned long depth,
                                        unsigned long minmatchlength,
                                        bool withsuftabcheck,
                                        bool countspms,
                                        bool outputspms)
{
  unsigned long idx, width, previous = 0;
  GtShortreadsortworkinfo *srsw;
  GtEncseqReader *esr, *esr1 = NULL, *esr2 = NULL;
  bool previousdefined = false;
  const uint16_t *lcptab_bucket;
  GtSpmsk_state *spmsk_state = NULL;
  unsigned long *suftab_bucket;

  if (withsuftabcheck)
  {
    esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
    esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  }
  if (outputspms || countspms)
  {
    spmsk_state = gt_spmsk_new(encseq,readmode,minmatchlength,
                               countspms,outputspms);
  }
  srsw = gt_shortreadsort_new(maxbucketsize,readmode,true);
  lcptab_bucket = gt_shortreadsort_lcpvalues(srsw);
  esr = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  suftab_bucket = gt_malloc(sizeof (*suftab_bucket) * maxbucketsize);
  for (idx = 0; idx <differentcodes; idx++)
  {
    gt_assert(countocc[idx+1]>countocc[idx]);
    width = countocc[idx+1] - countocc[idx];
    if (width >= 2UL)
    {
      gt_shortreadsort_array_sort(suftab_bucket,
                                  srsw,
                                  encseq,
                                  readmode,
                                  esr,
                                  spmsuftab,
                                  countocc[idx],
                                  width,
                                  depth);
      if (withsuftabcheck)
      {
        gt_firstcodes_checksuftab_bucket(encseq,
                                         readmode,
                                         esr1,
                                         esr2,
                                         previous,
                                         previousdefined,
                                         suftab_bucket,
                                         lcptab_bucket,
                                         width);
        previousdefined = true;
        previous = suftab_bucket[width-1];
      }
      if (outputspms || countspms)
      {
        int ret;

        ret = gt_spmsk_process(spmsk_state,
                               suftab_bucket,
                               lcptab_bucket,
                               width,
                               NULL);
        gt_assert(ret == 0);
      }
    }
  }
  gt_encseq_reader_delete(esr);
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  gt_shortreadsort_delete(srsw);
  gt_free(suftab_bucket);
  gt_spmsk_delete(spmsk_state);
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
                                                   unsigned int kmersize,
                                                   unsigned int minmatchlength,
                                                   bool withsuftabcheck,
                                                   bool countspms,
                                                   bool outputspms,
                                                   GtLogger *logger)
{
  GtTimer *timer = NULL;
  GtFirstcodesinfo fci;
  size_t sizeforcodestable, binsearchcache_size, suftab_size = 0;
  unsigned int numofchars;
  const unsigned int markprefixunits = 14U;
  const GtReadmode readmode = GT_READMODE_FORWARD;
  unsigned long totallength;
  unsigned long suftabentries;

  numofchars = gt_encseq_alphabetnumofchars(encseq);
  totallength = gt_encseq_total_length(encseq);
  fci.workspace = (size_t) gt_encseq_sizeofrep(encseq);
  fci.size_to_split = 0;
  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description("insert first codes into "
                                                   "array");
    gt_timer_start(timer);
  }
  fci.numofsequences = gt_encseq_num_of_sequences(encseq);
  gt_assert(numofchars == 4U);
  sizeforcodestable = sizeof (*fci.allfirstcodes) * fci.numofsequences;
  gt_logger_log(logger,"use array of size %lu",
                 (unsigned long) sizeforcodestable);
  fci.allfirstcodes = gt_malloc(sizeforcodestable);
  fci.differentcodes = 0;
  fci.countsequences = 0;
  fci.firstcodehits = 0;
  fci.firstcodeposhits = 0;
  fci.binsearchcache_unknownwidth = 0;
  fci.binsearchcache_unknowncount = 0;
  GT_INITARRAY(&fci.binsearchcache,GtIndexwithcode);
  getencseqkmers_twobitencoding(encseq,
                                readmode,
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
  fci.size_to_split = fci.markprefix.size;
  gt_marksubstring_init(&fci.marksuffix,numofchars,kmersize,
                        true,markprefixunits);
  fci.workspace += fci.marksuffix.size;
  fci.differentcodes = gt_remdups_in_sorted_array(&fci);
  gt_logger_log(logger,"number of different codes=%lu (%.4f) in %lu sequences",
                fci.differentcodes,
          (double) fci.differentcodes/fci.numofsequences,
          fci.countsequences);
  gt_logger_log(logger,"size of space to split %.1f megabytes",
                GT_MEGABYTES(fci.size_to_split));
  fci.binsearchcache_depth = (unsigned int) log10((double) fci.differentcodes);
  fci.flushcount = 0;
  fci.codebuffer_total = 0;
  binsearchcache_size = gt_firstcodes_halves(&fci,fci.binsearchcache_depth);
  fci.workspace += binsearchcache_size;
  gt_logger_log(logger,"binsearchcache_depth=%u => %lu bytes",
                       fci.binsearchcache_depth,
                       (unsigned long) binsearchcache_size);
  fci.codebuffer_allocated = fci.differentcodes/4;
  fci.codebuffer_nextfree = 0;
  fci.codebuffer = gt_malloc(sizeof (*fci.codebuffer)
                             * fci.codebuffer_allocated);
  fci.workspace += sizeof (*fci.codebuffer) * fci.codebuffer_allocated;
  gt_assert(fci.codebuffer_allocated > 0);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "accumulate counts",stdout);
  }
  fci.tempcodeforradixsort = gt_malloc((size_t) fci.codebuffer_allocated
                                       * sizeof (*fci.tempcodeforradixsort));
  fci.workspace += (size_t) fci.codebuffer_allocated
                        * sizeof (*fci.tempcodeforradixsort);
  getencseqkmers_twobitencoding(encseq,
                                readmode,
                                kmersize,
                                minmatchlength,
                                false,
                                gt_firstcodes_accumulatecounts,
                                &fci,
                                NULL,
                                NULL);
  gt_firstcodes_accumulatecounts_flush(&fci);
  gt_logger_log(logger,"codebuffer_total=%lu (%.2f)",
                fci.codebuffer_total,
                (double) fci.codebuffer_total/totallength);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "compute partial sums",stdout);
  }
  gt_logger_log(logger,"firstcodehits=%lu (%.2f)",fci.firstcodehits,
                                           (double) fci.firstcodehits/
                                           totallength);
  storefirstcodes_partialsum(&fci);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "insert suffixes",stdout);
  }
  gt_free(fci.tempcodeforradixsort);
  gt_free(fci.codebuffer);
  gt_assert(fci.codebuffer_nextfree == 0);
  fci.codebuffer_total = 0;
  fci.codebuffer_allocated /= 2UL;
  fci.codeposbuffer = gt_malloc(sizeof (*fci.codeposbuffer)
                                * fci.codebuffer_allocated);
  fci.tempcodeposforradixsort = gt_malloc(sizeof (*fci.tempcodeposforradixsort)
                                          * fci.codebuffer_allocated);
  suftabentries = fci.firstcodehits + fci.numofsequences;
  fci.spmsuftab = gt_spmsuftab_new(suftabentries,totallength);
  suftab_size = gt_spmsuftab_requiredspace(suftabentries,totallength);
  gt_logger_log(logger,"allocate %lu entries for suftab (%.2f megabytes)",
                suftabentries, GT_MEGABYTES(suftab_size));
  fci.flushcount = 0;
  getencseqkmers_twobitencoding(encseq,
                                readmode,
                                kmersize,
                                minmatchlength,
                                false,
                                gt_firstcodes_insertsuffixes,
                                &fci,
                                NULL,
                                NULL);
  gt_firstcodes_insertsuffixes_flush(&fci);
  gt_logger_log(logger,"firstcodeposhits=%lu",fci.firstcodeposhits);
  gt_assert(fci.firstcodeposhits == suftabentries);
  gt_logger_log(logger,"maxbucketsize=%lu",fci.maxbucketsize);
  GT_FREEARRAY(&fci.binsearchcache,GtIndexwithcode);
  gt_free(fci.codeposbuffer);
  gt_free(fci.tempcodeposforradixsort);
  gt_free(fci.markprefix.bits);
  gt_free(fci.marksuffix.bits);
  gt_free(fci.allfirstcodes);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "sort buckets of suffixes",stdout);
  }
  gt_firstcodes_sortremaining(encseq,
                              readmode,
                              fci.spmsuftab,
                              fci.maxbucketsize,
                              fci.countocc,
                              fci.differentcodes,
                              (unsigned long) kmersize,
                              (unsigned long) minmatchlength,
                              withsuftabcheck,
                              countspms,
                              outputspms);
  gt_free(fci.countocc);
  gt_spmsuftab_delete(fci.spmsuftab);
  gt_logger_log(logger,"workspace = %.2f",GT_MEGABYTES(fci.workspace));
  gt_logger_log(logger,"size to split = %.2f",GT_MEGABYTES(fci.size_to_split +
                                                           suftab_size));
  gt_logger_log(logger,"estimatedspace = %.2f",GT_MEGABYTES(fci.workspace +
                                                            fci.size_to_split +
                                                            suftab_size));
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
}
