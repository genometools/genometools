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
#include "firstcodes-tab.h"
#include "firstcodes.h"
#include "sfx-partssuf.h"

typedef struct
{
  unsigned long idx, code;
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

static unsigned long gt_kmercode_to_prefix_index(unsigned long code,
                                                 unsigned int shiftright)
{
  return code >> shiftright;
}

typedef struct
{
  unsigned long firstcodehits,
                firstcodeposhits,
                countsequences,
                numofsequences,
                codebuffer_allocated,
                codebuffer_nextfree,
                codebuffer_total,
                *codebuffer;
  GtUlongPair *codeposbuffer,
              *tempcodeposforradixsort;
  GtArrayGtIndexwithcode binsearchcache;
  unsigned int binsearchcache_depth,
               flushcount;
  Gtmarksubstring markprefix,
                  marksuffix;
  unsigned long *tempcodeforradixsort;
  GtSpmsuftab *spmsuftab;
  GtFirstcodestab tab;
} GtFirstcodesinfo;

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
    gt_marksubstring_mark(&fci->markprefix,fci->tab.allfirstcodes[0]);
    gt_marksubstring_mark(&fci->marksuffix,fci->tab.allfirstcodes[0]);
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
      gt_marksubstring_mark(&fci->markprefix,*readptr);
      gt_marksubstring_mark(&fci->marksuffix,*readptr);
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
                                           nextfreeGtIndexwithcode].idx = mid;
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

static unsigned long callfirstcodes_find = 0;

unsigned long gt_firstcodes_find(const GtFirstcodesinfo *fci,
                                 unsigned long code)
{
  unsigned long leftidx = ULONG_MAX, mididx, rightidx = ULONG_MAX;
  unsigned int depth;

  callfirstcodes_find++;
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
          gt_assert(leftic->idx != ULONG_MAX &&
                    rightic->idx != ULONG_MAX);
          if (leftic > fci->binsearchcache.spaceGtIndexwithcode)
          {
            leftidx = (leftic-1)->idx + 1;
          } else
          {
            leftidx = 0;
          }
          rightidx = rightic->idx - 1;
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
            gt_assert(leftic->idx != ULONG_MAX && rightic->idx != ULONG_MAX);
            leftidx = leftic->idx + 1;
            if (rightic < fci->binsearchcache.spaceGtIndexwithcode +
                          fci->binsearchcache.nextfreeGtIndexwithcode-1)
            {
              rightidx = (rightic+1)->idx - 1;
            } else
            {
              rightidx = fci->tab.differentcodes - 1;
            }
            break;
          }
        } else
        {
          return midic->idx;
        }
      }
    }
    gt_assert(leftidx != ULONG_MAX && rightidx != ULONG_MAX);
  } else
  {
    depth = 0;
    leftidx = 0;
    rightidx = fci->tab.differentcodes - 1;
  }
  while (leftidx <= rightidx)
  {
    mididx = leftidx + GT_DIV2(rightidx-leftidx);
    if (code < fci->tab.allfirstcodes[mididx])
    {
      rightidx = mididx - 1;
    } else
    {
      if (code > fci->tab.allfirstcodes[mididx])
      {
        leftidx = mididx + 1;
      } else
      {
        return mididx;
      }
    }
    depth++;
  }
  return ULONG_MAX;
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
                      *subjectstream_lst = fci->tab.allfirstcodes
                                           + fci->tab.differentcodes;

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
  unsigned long idx;
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
    idx = gt_firstcodes_find(fci,*vptr);
    if (idx != ULONG_MAX)
    {
      gt_assert(idx < fci->tab.differentcodes);
      fci->firstcodehits
        += gt_firstcodes_accumulatecounts_merge(fci,vptr,
                                                fci->tab.allfirstcodes + idx);
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
                      *subjectstream_lst = fci->tab.allfirstcodes +
                                           fci->tab.differentcodes;

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
  unsigned long idx;
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
       vptr < fci->codeposbuffer +
              fci->codebuffer_nextfree;
       vptr++)
  {
    idx = gt_firstcodes_find(fci,vptr->a);
    if (idx != ULONG_MAX)
    {
      gt_assert(idx < fci->tab.differentcodes);
      fci->firstcodeposhits
        += gt_firstcodes_insertsuffixes_merge(fci,vptr,
                                              fci->tab.allfirstcodes+idx);
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
  unsigned int numofchars;
  const unsigned int markprefixunits = 14U;
  const GtReadmode readmode = GT_READMODE_FORWARD;
  unsigned long maxbucketsize, totallength, suftabentries, size_to_split;
  GtSfxmappedrangelist *sfxmrlist = gt_Sfxmappedrangelist_new();
  size_t workspace;
  GtSuftabparts *suftabparts;

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
  gt_marksubstring_init(&fci.markprefix,numofchars,kmersize,false,
                        markprefixunits);
  gt_marksubstring_init(&fci.marksuffix,numofchars,kmersize,true,
                        markprefixunits);
  fci.tab.mappedmarkprefix = gt_Sfxmappedrange_new("markprefix",
                                                   fci.markprefix.entries,
                                                   GtSfxGtBitsequence,
                                                   gt_kmercode_to_prefix_index,
                                                   fci.markprefix.shiftright);
  gt_Sfxmappedrangelist_add(sfxmrlist,fci.tab.mappedmarkprefix);
  workspace += fci.marksuffix.size;
  fci.tab.differentcodes = gt_remdups_in_sorted_array(&fci);
  if (fci.tab.differentcodes > 0)
  {
    fci.tab.mappedallfirstcodes = gt_Sfxmappedrange_new("allfirstcodes",
                                                        fci.tab.differentcodes,
                                                        GtSfxunsignedlong,
                                                        NULL,0);
    gt_Sfxmappedrangelist_add(sfxmrlist,fci.tab.mappedallfirstcodes);
    fci.tab.mappedcountocc = gt_Sfxmappedrange_new("countocc",
                                                   fci.tab.differentcodes+1,
                                                   GtSfxunsignedlong,
                                                   NULL,0);
    gt_Sfxmappedrangelist_add(sfxmrlist,fci.tab.mappedcountocc);
  } else
  {
    fci.tab.mappedallfirstcodes = NULL;
    fci.tab.mappedcountocc = NULL;
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
  gt_assert(suftabparts != NULL);
  gt_suftabparts_showallrecords(suftabparts);
  gt_suftabparts_delete(suftabparts);
  gt_free(fci.tempcodeforradixsort);
  gt_free(fci.codebuffer);
  gt_assert(fci.codebuffer_nextfree == 0);
  fci.codebuffer_total = 0;
  fci.codebuffer_allocated /= 2UL;
  fci.codeposbuffer = gt_malloc(sizeof (*fci.codeposbuffer)
                                * fci.codebuffer_allocated);
  fci.tempcodeposforradixsort = gt_malloc(sizeof (*fci.tempcodeposforradixsort)
                                          * fci.codebuffer_allocated);
  fci.spmsuftab = gt_spmsuftab_new(suftabentries,totallength);
  suftab_size = gt_spmsuftab_requiredspace(suftabentries,totallength);
  gt_logger_log(logger,"allocate %lu entries for suftab (%.2f megabytes)",
                suftabentries, GT_MEGABYTES(suftab_size));
  fci.flushcount = 0;
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "insert suffixes",stdout);
  }
  /* Here the iteration over the parts starts */
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
  GT_FREEARRAY(&fci.binsearchcache,GtIndexwithcode);
  gt_free(fci.codeposbuffer);
  gt_free(fci.tempcodeposforradixsort);
  gt_free(fci.markprefix.bits);
  gt_free(fci.marksuffix.bits);
  gt_free(fci.tab.allfirstcodes);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "sort buckets of suffixes",stdout);
  }
  gt_firstcodes_sortremaining(encseq,
                              readmode,
                              fci.spmsuftab,
                              maxbucketsize,
                              fci.tab.countocc,
                              fci.tab.differentcodes,
                              (unsigned long) kmersize,
                              (unsigned long) minmatchlength,
                              withsuftabcheck,
                              countspms,
                              outputspms);
  /* Here the iteration over the parts ends */
  gt_free(fci.tab.countocc);
  gt_spmsuftab_delete(fci.spmsuftab);
  gt_Sfxmappedrangelist_delete(sfxmrlist);
  gt_Sfxmappedrange_delete(fci.tab.mappedmarkprefix,logger);
  gt_Sfxmappedrange_delete(fci.tab.mappedcountocc,logger);
  gt_Sfxmappedrange_delete(fci.tab.mappedallfirstcodes,logger);
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
  printf("callfirstcodes_find = %lu\n",callfirstcodes_find);
}
