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
#include "core/logger_api.h"
#include "core/spacepeak.h"
#include "core/spacecalc.h"
#include "core/fa.h"
#include "sfx-suffixer.h"
#include "sfx-shortreadsort.h"
#include "spmsuftab.h"
#include "esa-spmsk.h"
#include "firstcodes-tab.h"
#include "firstcodes.h"
#include "marksubstring.h"
#include "sfx-partssuf.h"

typedef struct
{
  unsigned long *ptr, code;
} GtIndexwithcode;

GT_DECLAREARRAYSTRUCT(GtIndexwithcode);

typedef struct
{
  unsigned long firstcodehits,
                firstcodeposhits,
                countsequences,
                numofsequences,
                codebuffer_allocated,
                codebuffer_nextfree,
                codebuffer_total,
                currentmincode,
                currentmaxcode,
                currentminindex,
                currentmaxindex,
                widthofpart,
                *codebuffer;
  GtUlongPair *codeposbuffer,
              *tempcodeposforradixsort;
  GtArrayGtIndexwithcode binsearchcache;
  unsigned int binsearchcache_depth,
               flushcount,
               shiftright2index;
  Gtmarksubstring *markprefix,
                  *marksuffix;
  unsigned long *tempcodeforradixsort;
  GtSpmsuftab *spmsuftab;
  GtSfxmappedrange *mappedcountocc,
                   *mappedallfirstcodes,
                   *mappedmarkprefix;
  GtFirstcodestab tab;
} GtFirstcodesinfo;

static unsigned long gt_kmercode_to_prefix_index(unsigned long idx,
                                                 const void *data)
{
  const GtFirstcodesinfo *fci = (const GtFirstcodesinfo *) data;

  return fci->tab.allfirstcodes[idx] >> fci->shiftright2index;
}

/* call the following function after computing the partial sums */

static void gt_storefirstcodes(void *processinfo,
                               GT_UNUSED bool firstinrange,
                               GT_UNUSED unsigned long pos,
                               GtCodetype code)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) processinfo;

  gt_assert(firstinrange);
  gt_assert(fci->tab.allfirstcodes != NULL &&
            fci->countsequences < fci->numofsequences);
  fci->tab.allfirstcodes[fci->countsequences++] = code;
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

    fci->tab.countocc = gt_calloc((size_t) (fci->numofsequences+1),
                              sizeof (*fci->tab.countocc));
    fci->tab.countocc[0] = 1UL;
    gt_marksubstring_mark(fci->markprefix,fci->tab.allfirstcodes[0]);
    gt_marksubstring_mark(fci->marksuffix,fci->tab.allfirstcodes[0]);
    for (storeptr = fci->tab.allfirstcodes, readptr = fci->tab.allfirstcodes+1;
         readptr < fci->tab.allfirstcodes + fci->numofsequences;
         readptr++)
    {
      if (*storeptr != *readptr)
      {
        storeptr++;

        *storeptr = *readptr;
      }
      fci->tab.countocc[(unsigned long) (storeptr - fci->tab.allfirstcodes)]++;
      gt_marksubstring_mark(fci->markprefix,*readptr);
      gt_marksubstring_mark(fci->marksuffix,*readptr);
    }
    numofdifferentcodes
      = (unsigned long) (storeptr - fci->tab.allfirstcodes + 1);
    if (numofdifferentcodes < fci->numofsequences)
    {
      /* reduce the memory requirement, as the duplicated elements are not
         needed */
      fci->tab.allfirstcodes = gt_realloc(fci->tab.allfirstcodes,
                                      sizeof (*fci->tab.allfirstcodes) *
                                      numofdifferentcodes);
      fci->tab.countocc = gt_realloc(fci->tab.countocc,
                                 sizeof (*fci->tab.countocc) *
                                 (numofdifferentcodes+1));
#ifdef SKDEBUG
      checkcodesorder(fci->tab.allfirstcodes,numofdifferentcodes,false);
#endif
    }
    gt_assert(fci->tab.countocc != NULL);
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
  }
  gt_assert(fci->binsearchcache.nextfreeGtIndexwithcode <
            fci->binsearchcache.allocatedGtIndexwithcode);
  fci->binsearchcache.spaceGtIndexwithcode[fci->binsearchcache.
                                           nextfreeGtIndexwithcode]
                                           .ptr = fci->tab.allfirstcodes + mid;
  fci->binsearchcache.spaceGtIndexwithcode[fci->binsearchcache.
                                           nextfreeGtIndexwithcode++]
                                           .code = fci->tab.allfirstcodes[mid];
  if (depth < maxdepth)
  {
    gt_firstcodes_halves_rek(fci,mid+1,right,depth+1,maxdepth);
  }
}

static size_t gt_firstcodes_halves(GtFirstcodesinfo *fci,
                                   unsigned int maxdepth)
{
  fci->binsearchcache.nextfreeGtIndexwithcode = 0;
  fci->binsearchcache.allocatedGtIndexwithcode = (1UL << (maxdepth+1)) - 1;
  if (fci->binsearchcache.allocatedGtIndexwithcode < fci->tab.differentcodes)
  {
    size_t allocbytes
      = sizeof (*fci->binsearchcache.spaceGtIndexwithcode)
                * fci->binsearchcache.allocatedGtIndexwithcode;
    fci->binsearchcache.spaceGtIndexwithcode = gt_malloc(allocbytes);
    gt_assert(fci->tab.differentcodes > 0);
    gt_firstcodes_halves_rek(fci,0,fci->tab.differentcodes - 1,0,maxdepth);
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
             fci->tab.allfirstcodes),
             fci->binsearchcache.spaceGtIndexwithcode[idx].code);
      }
    }
#endif
    return allocbytes;
  }
  return 0;
}

const unsigned long *gt_firstcodes_find(const GtFirstcodesinfo *fci,
                                        bool withcache,
                                        unsigned long leftbound,
                                        unsigned long rightbound,
                                        unsigned long code)
{
  const unsigned long *leftptr = NULL, *midptr, *rightptr = NULL;
  unsigned int depth;

  if (withcache && fci->binsearchcache.spaceGtIndexwithcode != NULL)
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
            leftptr = fci->tab.allfirstcodes;
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
              rightptr = fci->tab.allfirstcodes +
                         fci->tab.differentcodes - 1;
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
    leftptr = fci->tab.allfirstcodes + leftbound;
    rightptr = fci->tab.allfirstcodes + rightbound;
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
                                         + fci->codebuffer_nextfree - 1,
                      *subjectstream_lst = fci->tab.allfirstcodes
                                           + fci->tab.differentcodes - 1;

  while (query <= querystream_lst && subject <= subjectstream_lst)
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
        fci->tab.countocc[(unsigned long) (subject - fci->tab.allfirstcodes)]++;
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

  gt_assert(fci->tab.allfirstcodes != NULL);
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
    ptr = gt_firstcodes_find(fci,true,0,fci->tab.differentcodes-1,*vptr);
    if (ptr != NULL)
    {
      fci->firstcodehits += gt_firstcodes_accumulatecounts_merge(fci,vptr,ptr);
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
  GtCodetype tmpcode;

  if (!firstinrange &&
      GT_MARKSUBSTRING_CHECKMARK(fci->markprefix,code) &&
      GT_MARKSUBSTRING_CHECKMARK(fci->marksuffix,code))
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
                                       fci->codebuffer_nextfree - 1;
  const unsigned long *subject = subjectstream_fst,
                      *subjectstream_lst = fci->tab.allfirstcodes +
                                           fci->currentmaxindex;

  while (query <= querystream_lst && subject <= subjectstream_lst)
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
        idx = --fci->tab.countocc[(unsigned long)
                                  (subject - fci->tab.allfirstcodes)];
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

  gt_assert(fci->tab.allfirstcodes != NULL);
  gt_radixsort_GtUlongPair_linear(false,fci->codeposbuffer,
                                  fci->tempcodeposforradixsort,
                                  fci->codebuffer_nextfree);
#ifdef SKDEBUG
  checkcodeposorder(fci->codeposbuffer,fci->codebuffer_nextfree,true);
#endif
  fci->codebuffer_total += fci->codebuffer_nextfree;
  for (vptr = fci->codeposbuffer;
       vptr < fci->codeposbuffer + fci->codebuffer_nextfree;
       vptr++)
  {
    ptr = gt_firstcodes_find(fci,false,fci->currentminindex,
                                       fci->currentmaxindex,vptr->a);
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
  GtCodetype tmpcode;

  if (fci->currentmincode <= code &&
      code <= fci->currentmaxcode &&
      GT_MARKSUBSTRING_CHECKMARK(fci->markprefix,code) &&
      GT_MARKSUBSTRING_CHECKMARK(fci->marksuffix,code))
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

static unsigned long storefirstcodes_partialsum(GtFirstcodesinfo *fci)
{
  unsigned long idx, maxbucketsize;

  maxbucketsize = fci->tab.countocc[0];
  for (idx = 1UL; idx < fci->tab.differentcodes; idx++)
  {
    if (maxbucketsize < fci->tab.countocc[idx])
    {
      maxbucketsize = fci->tab.countocc[idx];
    }
    fci->tab.countocc[idx] += fci->tab.countocc[idx-1];
  }
  fci->tab.countocc[fci->tab.differentcodes]
    = fci->tab.countocc[fci->tab.differentcodes-1];
  return maxbucketsize;
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
                                        unsigned long minindex,
                                        unsigned long maxindex,
                                        unsigned long sumofwidth,
                                        unsigned long depth,
                                        GtSpmsk_state *spmsk_state,
                                        bool withsuftabcheck)
{
  unsigned long idx, width, previous = 0, sumwidth = 0;
  GtShortreadsortworkinfo *srsw;
  GtEncseqReader *esr, *esr1 = NULL, *esr2 = NULL;
  bool previousdefined = false;
  const uint16_t *lcptab_bucket;
  unsigned long *suftab_bucket;

  if (withsuftabcheck)
  {
    esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
    esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  }
  srsw = gt_shortreadsort_new(maxbucketsize,readmode,true);
  lcptab_bucket = gt_shortreadsort_lcpvalues(srsw);
  esr = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  suftab_bucket = gt_malloc(sizeof (*suftab_bucket) * maxbucketsize);
  for (idx = minindex; idx <=maxindex; idx++)
  {
    if (idx < maxindex)
    {
      gt_assert(countocc[idx+1]>countocc[idx]);
      width = countocc[idx+1] - countocc[idx];
    } else
    {
      gt_assert(sumofwidth > countocc[idx]);
      width = sumofwidth - countocc[idx];
    }
    sumwidth += width;
    gt_assert(sumwidth <= spmsuftab->numofentries);
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
      if (spmsk_state != NULL)
      {
        int ret;

        ret = gt_spmsk_process(spmsk_state,
                               suftab_bucket,
                               lcptab_bucket,
                               width,
                               NULL);
        gt_assert(ret == 0);
      }
    } else
    {
      gt_assert(width == 1UL);
    }
  }
  gt_encseq_reader_delete(esr);
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  gt_shortreadsort_delete(srsw);
  gt_free(suftab_bucket);
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
                                                   unsigned int parts,
                                                   unsigned int minmatchlength,
                                                   bool withsuftabcheck,
                                                   bool countspms,
                                                   bool outputspms,
                                                   GtLogger *logger)
{
  GtTimer *timer = NULL;
  GtFirstcodesinfo fci;
  size_t sizeforcodestable, binsearchcache_size, suftab_size = 0;
  unsigned int numofchars, part;
  const unsigned int markprefixunits = 14U;
  const GtReadmode readmode = GT_READMODE_FORWARD;
  unsigned long maxbucketsize, totallength, suftabentries, size_to_split,
                largest_width;
  GtSfxmappedrangelist *sfxmrlist = gt_Sfxmappedrangelist_new();
  size_t workspace;
  GtSuftabparts *suftabparts;
  GtSpmsk_state *spmsk_state = NULL;

  workspace = (size_t) gt_encseq_sizeofrep(encseq);
  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description("insert first codes into "
                                                   "array");
    gt_timer_start(timer);
  }
  fci.numofsequences = gt_encseq_num_of_sequences(encseq);
  sizeforcodestable = sizeof (*fci.tab.allfirstcodes) * fci.numofsequences;
  gt_logger_log(logger,"use array of size %lu",
                 (unsigned long) sizeforcodestable);
  fci.tab.allfirstcodes = gt_malloc(sizeforcodestable);
  fci.tab.differentcodes = 0;
  fci.countsequences = 0;
  fci.firstcodehits = 0;
  fci.firstcodeposhits = 0;
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
             fci.tab.allfirstcodes,fci.numofsequences);
  numofchars = gt_encseq_alphabetnumofchars(encseq);
  gt_assert(numofchars == 4U);
  fci.markprefix = gt_marksubstring_new(numofchars,kmersize,false,
                                        markprefixunits);
  fci.shiftright2index = gt_marksubstring_shiftright(fci.markprefix)
                         + GT_LOGWORDSIZE;
  fci.mappedmarkprefix
    = gt_Sfxmappedrange_new("markprefix",
                            gt_marksubstring_entries(fci.markprefix),
                            GtSfxGtBitsequence,
                            gt_kmercode_to_prefix_index,
                            &fci);
  gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedmarkprefix);
  fci.marksuffix = gt_marksubstring_new(numofchars,kmersize,true,
                                        markprefixunits);
  workspace += gt_marksubstring_size(fci.marksuffix);
  fci.tab.differentcodes = gt_remdups_in_sorted_array(&fci);
  if (fci.tab.differentcodes > 0)
  {
    fci.mappedallfirstcodes = gt_Sfxmappedrange_new("allfirstcodes",
                                                    fci.tab.differentcodes,
                                                    GtSfxunsignedlong,
                                                    NULL,NULL);
    gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedallfirstcodes);
    fci.mappedcountocc = gt_Sfxmappedrange_new("countocc",
                                               fci.tab.differentcodes+1,
                                               GtSfxunsignedlong,
                                               NULL,NULL);
    gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedcountocc);
  } else
  {
    fci.mappedallfirstcodes = NULL;
    fci.mappedcountocc = NULL;
  }
  size_to_split = gt_Sfxmappedrangelist_size_entire(sfxmrlist);
  gt_logger_log(logger,"number of different codes=%lu (%.4f) in %lu sequences",
                fci.tab.differentcodes,
                (double) fci.tab.differentcodes/fci.numofsequences,
                fci.countsequences);
  fci.binsearchcache_depth
    = (unsigned int) log10((double) fci.tab.differentcodes);
  fci.flushcount = 0;
  fci.codebuffer_total = 0;
  binsearchcache_size = gt_firstcodes_halves(&fci,fci.binsearchcache_depth);
  workspace += binsearchcache_size;
  gt_logger_log(logger,"binsearchcache_depth=%u => %lu bytes",
                       fci.binsearchcache_depth,
                       (unsigned long) binsearchcache_size);
  fci.codebuffer_allocated = fci.tab.differentcodes/4;
  fci.codebuffer_nextfree = 0;
  fci.codebuffer = gt_malloc(sizeof (*fci.codebuffer)
                             * fci.codebuffer_allocated);
  workspace += sizeof (*fci.codebuffer) * fci.codebuffer_allocated;
  gt_assert(fci.codebuffer_allocated > 0);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "accumulate counts",stdout);
  }
  fci.tempcodeforradixsort = gt_malloc((size_t) fci.codebuffer_allocated
                                       * sizeof (*fci.tempcodeforradixsort));
  workspace += (size_t) fci.codebuffer_allocated
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
  totallength = gt_encseq_total_length(encseq);
  gt_logger_log(logger,"codebuffer_total=%lu (%.2f)",
                fci.codebuffer_total,
                (double) fci.codebuffer_total/totallength);
  gt_logger_log(logger,"firstcodehits=%lu (%.2f)",fci.firstcodehits,
                          (double) fci.firstcodehits/totallength);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "compute partial sums",stdout);
  }
  maxbucketsize = storefirstcodes_partialsum(&fci);
  gt_logger_log(logger,"maxbucketsize=%lu",maxbucketsize);
  suftabentries = fci.firstcodehits + fci.numofsequences;
  gt_assert(parts > 0);
  suftabparts = gt_suftabparts_new(parts,
                                   NULL,
                                   &fci.tab,
                                   sfxmrlist,
                                   suftabentries,
                                   0,
                                   logger);
  gt_Sfxmappedrangelist_delete(sfxmrlist);
  sfxmrlist = NULL;
  gt_assert(suftabparts != NULL);
  /*
  gt_suftabparts_showallrecords(suftabparts);
  */
  gt_free(fci.tempcodeforradixsort);
  gt_free(fci.codebuffer);
  gt_assert(fci.codebuffer_nextfree == 0);
  if (gt_suftabparts_numofparts(suftabparts) > 1U)
  {
    gt_assert(fci.tab.allfirstcodes != NULL);
    gt_assert(fci.mappedallfirstcodes != NULL);
    gt_Sfxmappedrange_storetmp(fci.mappedallfirstcodes,
                               (void **) &fci.tab.allfirstcodes,
                               false,
                               logger);
    gt_assert(fci.tab.allfirstcodes == NULL);
    gt_assert(fci.tab.countocc != NULL);
    gt_assert(fci.mappedcountocc != NULL);
    gt_Sfxmappedrange_storetmp(fci.mappedcountocc,
                               (void **) &fci.tab.countocc,
                               true,
                               logger);
    gt_assert(fci.tab.countocc == NULL);

    gt_marksubstring_bits_null(fci.markprefix,false);
    gt_assert(fci.mappedmarkprefix != NULL);
    gt_Sfxmappedrange_storetmp(fci.mappedmarkprefix,
                               gt_marksubstring_bits_address(fci.markprefix),
                               false,
                               logger);
    gt_marksubstring_bits_null(fci.markprefix,true);
  } else
  {
    gt_Sfxmappedrange_delete(fci.mappedallfirstcodes,logger);
    fci.mappedallfirstcodes = NULL;
    gt_Sfxmappedrange_delete(fci.mappedcountocc,logger);
    fci.mappedcountocc = NULL;
    gt_Sfxmappedrange_delete(fci.mappedmarkprefix,logger);
    fci.mappedmarkprefix = NULL;
  }
  fci.codebuffer_total = 0;
  fci.codebuffer_allocated /= 2UL;
  fci.codeposbuffer = gt_malloc(sizeof (*fci.codeposbuffer)
                                * fci.codebuffer_allocated);
  fci.tempcodeposforradixsort = gt_malloc(sizeof (*fci.tempcodeposforradixsort)
                                          * fci.codebuffer_allocated);
  largest_width = gt_suftabparts_largest_width(suftabparts);
  fci.spmsuftab = gt_spmsuftab_new(largest_width,totallength,logger);
  suftab_size = gt_spmsuftab_requiredspace(largest_width,totallength);
  fci.flushcount = 0;
  if (outputspms || countspms)
  {
    spmsk_state = gt_spmsk_new(encseq,readmode,(unsigned long ) minmatchlength,
                               countspms,outputspms);
  }
  for (part = 0; part < gt_suftabparts_numofparts(suftabparts); part++)
  {
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "to insert suffixes into buckets",stdout);
    }
    gt_logger_log(logger,"compute part %u",part);
    fci.currentminindex = gt_suftabparts_minindex(part,suftabparts);
    fci.currentmaxindex = gt_suftabparts_maxindex(part,suftabparts);
    fci.widthofpart = gt_suftabparts_widthofpart(part,suftabparts);
    if (fci.mappedallfirstcodes != NULL)
    {
      fci.tab.allfirstcodes
        = (unsigned long *)
          gt_Sfxmappedrange_map(fci.mappedallfirstcodes,
                                part,
                                fci.currentminindex,
                                fci.currentmaxindex,
                                logger);
    }
    if (fci.mappedcountocc != NULL)
    {
      fci.tab.countocc
        = (unsigned long *)
          gt_Sfxmappedrange_map(fci.mappedcountocc,
                                part,
                                fci.currentminindex,
                                fci.currentmaxindex,
                                logger);
    }
    if (fci.mappedmarkprefix != NULL)
    {
      gt_marksubstring_bits_map(fci.markprefix,
                                (GtBitsequence *)
                                gt_Sfxmappedrange_map(fci.mappedmarkprefix,
                                                      part,
                                                      fci.currentminindex,
                                                      fci.currentmaxindex,
                                                      logger));
    }
    fci.currentmincode = gt_suftabparts_mincode(part,suftabparts);
    fci.currentmaxcode = gt_suftabparts_maxcode(part,suftabparts);
    gt_spmsuftab_partoffset(fci.spmsuftab,
                            gt_suftabparts_offset(part,suftabparts));
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
    if (part == gt_suftabparts_numofparts(suftabparts) - 1)
    {
      GT_FREEARRAY(&fci.binsearchcache,GtIndexwithcode);
      gt_free(fci.codeposbuffer);
      gt_free(fci.tempcodeposforradixsort);
      gt_marksubstring_delete(fci.markprefix,
                              fci.mappedmarkprefix == NULL ? true : false);
      gt_marksubstring_delete(fci.marksuffix,true);
      if (fci.mappedallfirstcodes == NULL)
      {
        gt_free(fci.tab.allfirstcodes);
      }
    }
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "sort buckets of suffixes",stdout);
    }
    gt_firstcodes_sortremaining(encseq,
                                readmode,
                                fci.spmsuftab,
                                maxbucketsize,
                                fci.tab.countocc,
                                gt_suftabparts_minindex(part,suftabparts),
                                gt_suftabparts_maxindex(part,suftabparts),
                                gt_suftabparts_sumofwidth(part,suftabparts),
                                (unsigned long) kmersize,
                                spmsk_state,
                                withsuftabcheck);
  }
  gt_spmsk_delete(spmsk_state);
  gt_suftabparts_delete(suftabparts);
  gt_logger_log(logger,"firstcodeposhits=%lu",fci.firstcodeposhits);
  gt_assert(fci.firstcodeposhits == suftabentries);
  if (fci.mappedcountocc == NULL)
  {
    gt_free(fci.tab.countocc);
  }
  gt_spmsuftab_delete(fci.spmsuftab);
  gt_Sfxmappedrange_delete(fci.mappedmarkprefix,logger);
  gt_Sfxmappedrange_delete(fci.mappedcountocc,logger);
  gt_Sfxmappedrange_delete(fci.mappedallfirstcodes,logger);
  gt_logger_log(logger,"workspace = %.2f",GT_MEGABYTES(workspace));
  gt_logger_log(logger,"size to split = %.2f",
                GT_MEGABYTES(size_to_split + suftab_size));
  gt_logger_log(logger,"estimatedspace = %.2f",
                GT_MEGABYTES(workspace +
                             size_to_split +
                             suftab_size));
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  /*printf("callfirstcodes_find = %lu\n",callfirstcodes_find);*/
}
