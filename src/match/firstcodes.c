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
#include "core/error_api.h"
#include "sfx-suffixer.h"
#include "sfx-shortreadsort.h"
#include "spmsuftab.h"
#include "firstcodes-tab.h"
#include "firstcodes-buf.h"
#include "firstcodes.h"
#include "marksubstring.h"
#include "sfx-partssuf.h"
#include "seqnumrelpos.h"

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
                codebuffer_total,
                currentminindex,
                currentmaxindex,
                widthofpart;
  GtUlongPair *tempcodeposforradixsort;
  GtArrayGtIndexwithcode binsearchcache;
  unsigned int binsearchcache_depth,
               flushcount,
               shiftright2index;
  unsigned long *tempcodeforradixsort;
  GtSpmsuftab *spmsuftab;
  GtSfxmappedrange *mappedcountocc,
                   *mappedoverflow,
                   *mappedallfirstcodes,
                   *mappedmarkprefix;
  GtCodeposbuffer buf;
  GtFirstcodestab tab;
} GtFirstcodesinfo;

static unsigned long gt_kmercode_to_prefix_index(unsigned long idx,
                                                 const void *data)
{
  const GtFirstcodesinfo *fci = (const GtFirstcodesinfo *) data;

  gt_assert(idx < fci->tab.differentcodes);
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

    gt_firstcodes_countocc_new(&fci->tab,fci->numofsequences);
    gt_firstcodes_countocc_increment(&fci->tab,0);
    gt_marksubstring_mark(fci->buf.markprefix,fci->tab.allfirstcodes[0]);
    gt_marksubstring_mark(fci->buf.marksuffix,fci->tab.allfirstcodes[0]);
    for (storeptr = fci->tab.allfirstcodes, readptr = fci->tab.allfirstcodes+1;
         readptr < fci->tab.allfirstcodes + fci->numofsequences;
         readptr++)
    {
      if (*storeptr != *readptr)
      {
        storeptr++;

        *storeptr = *readptr;
      }
      gt_firstcodes_countocc_increment(&fci->tab,(unsigned long)
                                       (storeptr - fci->tab.allfirstcodes));
      gt_marksubstring_mark(fci->buf.markprefix,*readptr);
      gt_marksubstring_mark(fci->buf.marksuffix,*readptr);
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
      gt_firstcodes_countocc_resize(&fci->tab,numofdifferentcodes);
#ifdef SKDEBUG
      checkcodesorder(fci->tab.allfirstcodes,numofdifferentcodes,false);
#endif
    }
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
                      *querystream_lst = fci->buf.spaceGtUlong
                                         + fci->buf.nextfree - 1,
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
        gt_firstcodes_countocc_increment(&fci->tab,(unsigned long)
                                         (subject - fci->tab.allfirstcodes));
        query++;
        found++;
      }
    }
  }
  return found;
}

static void gt_firstcodes_accumulatecounts_flush(void *data)
{
  const unsigned long *ptr;
  unsigned long *vptr;
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) data;

  gt_assert(fci->tab.allfirstcodes != NULL);
  gt_radixsort_GtUlong_linear(false,fci->buf.spaceGtUlong,
                              fci->tempcodeforradixsort,
                              fci->buf.nextfree);
#ifdef SKDEBUG
  checkcodesorder(fci->buf.spaceGtUlong,fci->buf.nextfree,true);
#endif
  fci->codebuffer_total += fci->buf.nextfree;
  for (vptr = fci->buf.spaceGtUlong;
       vptr < fci->buf.spaceGtUlong + fci->buf.nextfree;
       vptr++)
  {
    ptr = gt_firstcodes_find(fci,true,0,fci->tab.differentcodes-1,*vptr);
    if (ptr != NULL)
    {
      fci->firstcodehits += gt_firstcodes_accumulatecounts_merge(fci,vptr,ptr);
      break;
    }
  }
  fci->flushcount++;
  fci->buf.nextfree = 0;
}

static unsigned long gt_firstcodes_insertsuffixes_merge(
                                        GtFirstcodesinfo *fci,
                                        const GtUlongPair *querystream_fst,
                                        const unsigned long *subjectstream_fst)
{
  unsigned long found = 0, idx;
  const GtUlongPair *query = querystream_fst,
                    *querystream_lst = fci->buf.spaceGtUlongPair +
                                       fci->buf.nextfree - 1;
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
        idx = (unsigned long) (subject - fci->tab.allfirstcodes);
        idx = gt_firstcodes_insertionindex(&fci->tab,idx);
        gt_assert(idx < fci->firstcodehits + fci->numofsequences);
        gt_spmsuftab_set(fci->spmsuftab,idx,
                         gt_spmsuftab_usebitsforpositions(fci->spmsuftab)
                           ? gt_seqnumrelpos_decode_pos(fci->buf.snrp,query->b)
                           : query->b);
        query++;
        found++;
      }
    }
  }
  return found;
}

static void gt_firstcodes_insertsuffixes_flush(void *data)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) data;
  const unsigned long *ptr;
  GtUlongPair *vptr;

  gt_assert(fci->tab.allfirstcodes != NULL);
  gt_radixsort_GtUlongPair_linear(false,fci->buf.spaceGtUlongPair,
                                  fci->tempcodeposforradixsort,
                                  fci->buf.nextfree);
#ifdef SKDEBUG
  checkcodeposorder(fci->buf.spaceGtUlongPair,fci->buf.nextfree,true);
#endif
  fci->codebuffer_total += fci->buf.nextfree;
  for (vptr = fci->buf.spaceGtUlongPair;
       vptr < fci->buf.spaceGtUlongPair + fci->buf.nextfree;
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
  fci->flushcount++;
  fci->buf.nextfree = 0;
}

static void gt_firstcodes_checksuftab_bucket(const GtEncseq *encseq,
                                             GtReadmode readmode,
                                             GtEncseqReader *esr1,
                                             GtEncseqReader *esr2,
                                             unsigned long previoussuffix,
                                             bool previousdefined,
                                             const unsigned long
                                               *seqnum_relpos_bucket,
                                             const GtSeqnumrelpos *snrp,
                                             GT_UNUSED const uint16_t
                                               *lcptab_bucket,
                                             unsigned long numberofsuffixes)
{
  unsigned long idx, current, maxlcp,
                totallength = gt_encseq_total_length(encseq);
  int cmp;
  const bool specialsareequal = false, specialsareequalatdepth0 = false;
  const unsigned long depth = 0;

  gt_assert(!previousdefined || previoussuffix < totallength);
  for (idx = 0; idx < numberofsuffixes; idx++)
  {
    current = gt_seqnumrelpos_decode_pos(snrp,seqnum_relpos_bucket[idx]);
    if (previousdefined && idx < totallength)
    {
      gt_assert(current < totallength);
      cmp = gt_encseq_check_comparetwosuffixes(encseq,
                                               readmode,
                                               &maxlcp,
                                               specialsareequal,
                                               specialsareequalatdepth0,
                                               depth,
                                               previoussuffix,
                                               current,
                                               esr1,
                                               esr2);
      gt_assert(cmp <= 0);
      gt_assert(idx == 0 || maxlcp == (unsigned long) lcptab_bucket[idx]);
    }
    previoussuffix = current;
    previousdefined = true;
  }
}

static int gt_firstcodes_sortremaining(const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       GtShortreadsortworkinfo *srsw,
                                       unsigned long *seqnum_relpos_bucket,
                                       GtSpmsuftab *spmsuftab,
                                       const GtSeqnumrelpos *snrp,
                                       const GtFirstcodestab *fct,
                                       unsigned long minindex,
                                       unsigned long maxindex,
                                       unsigned long sumofwidth,
                                       unsigned long depth,
                                       GtFirstcodesintervalprocess itvprocess,
                                       void *itvprocessdata,
                                       bool withsuftabcheck,
                                       GtError *err)
{
  unsigned long current, next, idx, width, previoussuffix = 0, sumwidth = 0;
  GtEncseqReader *esr1 = NULL, *esr2 = NULL;
  bool previousdefined = false;
  const uint16_t *lcptab_bucket;
  bool haserr = false;

  if (withsuftabcheck)
  {
    esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
    esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  }
  lcptab_bucket = gt_shortreadsort_lcpvalues(srsw);
  current = gt_firstcodes_get_leftborder(fct,minindex);
  for (idx = minindex; idx <=maxindex; idx++)
  {
    if (idx < maxindex)
    {
      next = gt_firstcodes_get_leftborder(fct,idx+1);
      gt_assert(next > current);
      width = next - current;
    } else
    {
      gt_assert(sumofwidth > current);
      width = sumofwidth - current;
    }
    sumwidth += width;
    gt_assert(sumwidth <= spmsuftab->numofentries);
    if (width >= 2UL)
    {
      gt_shortreadsort_array_sort(seqnum_relpos_bucket,
                                  snrp,
                                  srsw,
                                  encseq,
                                  spmsuftab,
                                  current,
                                  width,
                                  depth);
      if (withsuftabcheck)
      {
        gt_firstcodes_checksuftab_bucket(encseq,
                                         readmode,
                                         esr1,
                                         esr2,
                                         previoussuffix,
                                         previousdefined,
                                         seqnum_relpos_bucket,
                                         snrp,
                                         lcptab_bucket,
                                         width);
        previousdefined = true;
        previoussuffix
          = gt_seqnumrelpos_decode_pos(snrp,seqnum_relpos_bucket[width-1]);
      }
      if (itvprocess != NULL)
      {
        if (itvprocess(itvprocessdata,seqnum_relpos_bucket,snrp,
                        lcptab_bucket, width, err) != 0)
        {
          haserr = true;
          break;
        }
      }
    } else
    {
      gt_assert(width == 1UL);
    }
    current = next;
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  return haserr ? -1 : 0;
}

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME)                        firstcodes_##NAME
#define firstcodes_ARRAY_GET(ARR,RELIDX)       ARR[RELIDX]
#define firstcodes_ARRAY_SET(ARR,RELIDX,VALUE) ARR[RELIDX] = VALUE

typedef unsigned long QSORTNAME(Sorttype);

#include "match/qsort-direct.gen"

typedef struct
{
  size_t workspace, splitspace;
} GtFirstcodesspacelog;

static void gt_firstcodes_update_workspace(int line,
                                           bool add,
                                           const char *kind,
                                           bool addtowork,
                                           GtFirstcodesspacelog *fcsl,
                                           size_t workspace)
{
  if (addtowork)
  {
    if (add)
    {
      fcsl->workspace += workspace;
    } else
    {
      gt_assert(fcsl->workspace >= workspace);
      fcsl->workspace -= workspace;
    }
  } else
  {
    if (add)
    {
      fcsl->splitspace += workspace;
    } else
    {
      gt_assert(fcsl->workspace >= workspace);
      fcsl->splitspace -= workspace;
    }
  }
  gt_log_log("line %d: %s %.2f MB for %s to %space; work=%.2f, "
             "split=%.2f,all=%.2f MB",
             line,
             add ? "add" : "delete",
             GT_MEGABYTES(workspace),
             kind,
             addtowork ? "work" : "split",
             GT_MEGABYTES(fcsl->workspace),
             GT_MEGABYTES(fcsl->splitspace),
             GT_MEGABYTES(fcsl->workspace+fcsl->splitspace));
  /*
  gt_logger_log(logger,"current space peak %.2f",
                GT_MEGABYTES(gt_spacepeak_get_space_peak()));
  */
}

static void gt_firstcode_delete_before_end(GtFirstcodesinfo *fci)
{
  GT_FREEARRAY(&fci->binsearchcache,GtIndexwithcode);
  gt_free(fci->buf.spaceGtUlongPair);
  gt_free(fci->tempcodeposforradixsort);
  gt_Sfxmappedrange_delete(fci->mappedmarkprefix);
  gt_marksubstring_delete(fci->buf.markprefix,true);
  fci->buf.markprefix = NULL;
  gt_marksubstring_delete(fci->buf.marksuffix,true);
  if (fci->mappedallfirstcodes == NULL)
  {
    gt_free(fci->tab.allfirstcodes);
  }
}

int storefirstcodes_getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                                   unsigned int kmersize,
                                                   unsigned int numofparts,
                                                   unsigned long maximumspace,
                                                   unsigned int minmatchlength,
                                                   bool withsuftabcheck,
                                                   bool onlyaccumulation,
                                                   GtFirstcodesintervalprocess
                                                     itvprocess,
                                                   void *itvprocessdata,
                                                   GtLogger *logger,
                                                   GtError *err)
{
  GtTimer *timer = NULL;
  GtFirstcodesinfo fci;
  size_t sizeforcodestable, binsearchcache_size, suftab_size = 0;
  unsigned int numofchars, part, bitsforrelpos, bitsforseqnum;
  const unsigned int markprefixunits = 14U;
  const unsigned int marksuffixunits = 14U;
  const GtReadmode readmode = GT_READMODE_FORWARD;
  unsigned long maxbucketsize, maxseqlength, numofdbsequences, maxrelpos,
                totallength, suftabentries = 0, largest_width;
  GtSfxmappedrangelist *sfxmrlist;
  GtSuftabparts *suftabparts = NULL, *suftabparts_sample;
  GtFirstcodesspacelog fcsl;
  GtShortreadsortworkinfo *srsw = NULL;
  unsigned long *seqnum_relpos_bucket = NULL;
  bool haserr = false;

  maxseqlength = gt_encseq_max_seq_length(encseq);
  totallength = gt_encseq_total_length(encseq);
  if (maxseqlength > (unsigned long) minmatchlength)
  {
    maxrelpos = maxseqlength - (unsigned long) minmatchlength;
  } else
  {
    maxrelpos = 0;
  }
  bitsforrelpos = gt_determinebitspervalue((uint64_t) maxrelpos);
  fci.buf.snrp = gt_seqnumrelpos_new(bitsforrelpos,encseq);
  numofdbsequences = gt_encseq_num_of_sequences(encseq);
  gt_assert(numofdbsequences > 0);
  bitsforseqnum = gt_determinebitspervalue((uint64_t) (numofdbsequences - 1));
  if (bitsforseqnum + bitsforrelpos > (unsigned int) GT_INTWORDSIZE)
  {
    gt_error_set(err,"cannot process encoded sequences with %lu sequences "
                     "of length up to %lu (%u+%u bits)",
                     numofdbsequences,maxseqlength,bitsforseqnum,bitsforrelpos);
    return -1;
  }
  sfxmrlist = gt_Sfxmappedrangelist_new();
  fcsl.workspace = 0;
  fcsl.splitspace = 0;
  fci.spmsuftab = NULL;
  fci.buf.spaceGtUlongPair = NULL;
  fci.buf.spaceGtUlong = NULL;
  fci.tempcodeposforradixsort = NULL;
  fci.tempcodeforradixsort = NULL;
  fci.mappedoverflow = NULL;
  gt_firstcodes_update_workspace(__LINE__,true,"encseq",true,&fcsl,
                                 (size_t) gt_encseq_sizeofrep(encseq));
  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description("insert first codes into "
                                                   "array");
    gt_timer_start(timer);
  }
  fci.numofsequences = gt_encseq_num_of_sequences(encseq);
  sizeforcodestable = sizeof (*fci.tab.allfirstcodes) * fci.numofsequences;
  gt_logger_log(logger,"store %lu prefix codes",fci.numofsequences);
  fci.tab.allfirstcodes = gt_malloc(sizeforcodestable);
  gt_firstcodes_countocc_setnull(&fci.tab);
  fci.tab.overflow_leftborder = NULL;
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
    gt_timer_show_progress(timer, "to sort the codes",stdout);
  }
  QSORTNAME(gt_direct_qsort)
             (6UL, false,
             fci.tab.allfirstcodes,fci.numofsequences);
  numofchars = gt_encseq_alphabetnumofchars(encseq);
  gt_assert(numofchars == 4U);
  fci.buf.markprefix = gt_marksubstring_new(numofchars,kmersize,false,
                                            markprefixunits);
  fci.shiftright2index = gt_marksubstring_shiftright(fci.buf.markprefix)
                         + GT_LOGWORDSIZE;
  fci.mappedmarkprefix
    = gt_Sfxmappedrange_new("markprefix",
                            gt_marksubstring_entries(fci.buf.markprefix),
                            GtSfxGtBitsequence,
                            gt_kmercode_to_prefix_index,
                            &fci);
  gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedmarkprefix);
  fci.buf.marksuffix = gt_marksubstring_new(numofchars,kmersize,true,
                                            marksuffixunits);
  gt_firstcodes_update_workspace(__LINE__,true,"marksuffix",true,&fcsl,
                                 (size_t)
                                 gt_marksubstring_size(fci.buf.marksuffix));
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
                                               GtSfxuint32_t,
                                               NULL,NULL);
    gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedcountocc);
  } else
  {
    fci.mappedallfirstcodes = NULL;
    fci.mappedcountocc = NULL;
  }
  gt_firstcodes_update_workspace(__LINE__,true,"splitable arrays",false,&fcsl,
                                 (size_t)
                                 gt_Sfxmappedrangelist_size_entire(sfxmrlist));
  gt_logger_log(logger,"number of different codes=%lu (%.4f) in %lu sequences",
                fci.tab.differentcodes,
                (double) fci.tab.differentcodes/fci.numofsequences,
                fci.countsequences);
  fci.binsearchcache_depth
    = (unsigned int) log10((double) fci.tab.differentcodes);
  fci.flushcount = 0;
  fci.codebuffer_total = 0;
  binsearchcache_size = gt_firstcodes_halves(&fci,fci.binsearchcache_depth);
  gt_firstcodes_update_workspace(__LINE__,true,"binsearchcache",
                                 true,&fcsl,binsearchcache_size);
  gt_log_log("binsearchcache_depth=%u => %lu bytes",
             fci.binsearchcache_depth,(unsigned long) binsearchcache_size);
  if (maximumspace > 0)
  {
    if ((unsigned long) (fcsl.workspace + fcsl.splitspace) >= maximumspace)
    {
      gt_error_set(err,"already used %.2f MB of memory, cannot compute "
                       "index in at most %.2f MB",
                       GT_MEGABYTES(fcsl.workspace + fcsl.splitspace),
                       GT_MEGABYTES(maximumspace));
      haserr = true;
    } else
    {
      size_t remainspace = (size_t) maximumspace -
                           (fcsl.workspace + fcsl.splitspace);

      fci.buf.allocated = (unsigned long)
                          remainspace / (sizeof (*fci.buf.spaceGtUlong) +
                                         sizeof (*fci.tempcodeforradixsort));
      if (fci.buf.allocated < fci.tab.differentcodes/16UL)
      {
        fci.buf.allocated = fci.tab.differentcodes/16UL;
      }
    }
  } else
  {
    fci.buf.allocated = fci.tab.differentcodes/5;
  }
  if (fci.buf.allocated < 16UL)
  {
    fci.buf.allocated = 16UL;
  }
  fci.buf.nextfree = 0;
  if (!haserr)
  {
    fci.buf.spaceGtUlong = gt_malloc(sizeof (*fci.buf.spaceGtUlong)
                                     * fci.buf.allocated);
    gt_firstcodes_update_workspace(__LINE__,true,"position buffer",true,&fcsl,
                                   sizeof (*fci.buf.spaceGtUlong)
                                   * fci.buf.allocated);
    gt_assert(fci.buf.allocated > 0);
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "to accumulate counts",stdout);
    }
    fci.tempcodeforradixsort = gt_malloc((size_t) fci.buf.allocated
                                         * sizeof (*fci.tempcodeforradixsort));
    gt_firstcodes_update_workspace(__LINE__,true,"tempcodeforradixsort",
                                   true,&fcsl,
                                   (size_t) fci.buf.allocated
                                   * sizeof (*fci.tempcodeforradixsort));

    fci.buf.fciptr = &fci; /* as we need to give fci to the flush function */
    fci.buf.flush_function = gt_firstcodes_accumulatecounts_flush;
    gt_firstcodes_accumulatecounts_getencseqkmers_twobitencoding(
                                  encseq,
                                  readmode,
                                  kmersize,
                                  minmatchlength,
                                  &fci.buf,
                                  NULL);
    gt_firstcodes_accumulatecounts_flush(&fci);
    gt_logger_log(logger,"codebuffer_total=%lu (%.3f%% of all suffixes)",
                  fci.codebuffer_total,
                  100.0 * (double) fci.codebuffer_total/totallength);
    gt_logger_log(logger,"firstcodehits=%lu (%.3f%% of all suffixes), "
                         "%u rounds (avg length %lu)",
                           fci.firstcodehits,
                           100.0 * (double) fci.firstcodehits/totallength,
                           fci.flushcount,
                           fci.firstcodehits/fci.flushcount);
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "to compute partial sums",stdout);
    }
    maxbucketsize = gt_firstcodes_partialsums(&fci.tab);
    gt_logger_log(logger,"maxbucketsize=%lu",maxbucketsize);
    if (fci.tab.overflow_index > 0)
    {
      unsigned long overflowcells = fci.tab.differentcodes -
                                    fci.tab.overflow_index + 1;
      fci.mappedoverflow = gt_Sfxmappedrange_new("overflow_leftborder",
                                                 overflowcells,
                                                 GtSfxunsignedlong,NULL,NULL);
      gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedoverflow);
    }
    suftabentries = fci.firstcodehits + fci.numofsequences;
    if (maximumspace > 0)
    {
      int retval;

      gt_assert(numofparts == 1U);
      retval = gt_suftabparts_fit_memlimit(fcsl.workspace + fcsl.splitspace,
                                           maximumspace,
                                           NULL,
                                           &fci.tab,
                                           sfxmrlist,
                                           totallength,
                                           bitsforseqnum + bitsforrelpos,
                                           0, /* special characters not used */
                                           suftabentries,
                                           false, /* suftabuint not used */
                                           err);
      if (retval < 0)
      {
        haserr = true;
      } else
      {
        gt_assert(retval > 0);
        numofparts = (unsigned int) retval;
        gt_logger_log(logger, "derived parts=%u",numofparts);
      }
    }
  }
  if (!haserr)
  {
    gt_assert(numofparts > 0);
    suftabparts = gt_suftabparts_new(numofparts,
                                     NULL,
                                     &fci.tab,
                                     sfxmrlist,
                                     suftabentries,
                                     0,
                                     logger);
    gt_suftabparts_showallrecords(suftabparts,true);
    fci.tab.usesample = true;
    suftabparts_sample = gt_suftabparts_new(numofparts,
                                      NULL,
                                      &fci.tab,
                                      sfxmrlist,
                                      suftabentries,
                                      0,
                                      logger);
    gt_suftabparts_showallrecords(suftabparts_sample,true);
    gt_assert(fci.tab.allfirstcodes[fci.tab.differentcodes - 1] ==
              gt_suftabparts_maxindex_last(suftabparts_sample));
    gt_suftabparts_delete(suftabparts_sample);
    fci.tab.usesample = false;
    gt_Sfxmappedrangelist_delete(sfxmrlist);
    sfxmrlist = NULL;
    gt_assert(suftabparts != NULL);
    /*
    gt_suftabparts_showallrecords(suftabparts);
    */
    gt_free(fci.tempcodeforradixsort);
    gt_firstcodes_update_workspace(__LINE__,false,"tempcodeforradixsort",
                                   true,&fcsl,
                                   (size_t) fci.buf.allocated
                                   * sizeof (*fci.tempcodeforradixsort));
    gt_free(fci.buf.spaceGtUlong);
    gt_firstcodes_update_workspace(__LINE__,false,"buf.spaceGtUlong",
                                   true,&fcsl,
                                   (size_t) fci.buf.allocated
                                   * sizeof (*fci.buf.spaceGtUlong));
    gt_assert(fci.buf.nextfree == 0);
    if (gt_suftabparts_numofparts(suftabparts) > 1U)
    {
      gt_assert(fci.tab.allfirstcodes != NULL);
      gt_assert(fci.mappedallfirstcodes != NULL);
      gt_Sfxmappedrange_storetmp(fci.mappedallfirstcodes,
                                 (void **) &fci.tab.allfirstcodes,
                                 false);
      gt_assert(fci.tab.allfirstcodes == NULL);
      gt_assert(fci.mappedcountocc != NULL);
      gt_Sfxmappedrange_storetmp(fci.mappedcountocc,
                                 gt_firstcodes_countocc_address(&fci.tab),
                                 true);
      if (fci.tab.overflow_index > 0)
      {
        gt_assert(fci.mappedoverflow != NULL);
        gt_Sfxmappedrange_storetmp(fci.mappedoverflow,
                                   (void **) &fci.tab.overflow_leftborder,
                                   true);
      }
      gt_marksubstring_bits_null(fci.buf.markprefix,false);
      gt_assert(fci.mappedmarkprefix != NULL);
      gt_Sfxmappedrange_storetmp(fci.mappedmarkprefix,
                                 gt_marksubstring_bits_address(
                                    fci.buf.markprefix),
                                 false);
      gt_firstcodes_update_workspace(__LINE__,false,"splitable arrays",
                                     false,&fcsl,
                                     (size_t)
                                     fcsl.splitspace);
      gt_marksubstring_bits_null(fci.buf.markprefix,true);
    } else
    {
      gt_Sfxmappedrange_delete(fci.mappedallfirstcodes);
      fci.mappedallfirstcodes = NULL;
      gt_Sfxmappedrange_delete(fci.mappedcountocc);
      fci.mappedcountocc = NULL;
      gt_Sfxmappedrange_delete(fci.mappedoverflow);
      fci.mappedoverflow = NULL;
      gt_Sfxmappedrange_delete(fci.mappedmarkprefix);
      fci.mappedmarkprefix = NULL;
    }
    fci.codebuffer_total = 0;
    largest_width = gt_suftabparts_largest_width(suftabparts);
    fci.spmsuftab = gt_spmsuftab_new(largest_width,totallength,
                                     bitsforseqnum + bitsforrelpos,
                                     logger);
    suftab_size = gt_spmsuftab_requiredspace(largest_width,totallength,
                                             bitsforseqnum + bitsforrelpos);
    fci.buf.flush_function = gt_firstcodes_insertsuffixes_flush;
    srsw = gt_shortreadsort_new(maxbucketsize,readmode,true);
    gt_firstcodes_update_workspace(__LINE__,true,"shortreadsort",
                                   true,&fcsl,
                                   (size_t)
                                   gt_shortreadsort_size(true,maxbucketsize));
    seqnum_relpos_bucket
      = gt_malloc(sizeof (*seqnum_relpos_bucket) * maxbucketsize);
    gt_firstcodes_update_workspace(__LINE__,true,"seqnum_relpos_bucket",
                                   true,&fcsl,
                                   (size_t)
                                   sizeof (*seqnum_relpos_bucket) *
                                       maxbucketsize);
    if (maximumspace > 0)
    {
      size_t used = fcsl.workspace + suftab_size +
                    gt_suftabparts_largestsizemappedpartwise(suftabparts);

      gt_assert(maximumspace > (unsigned long) used);
      fci.buf.allocated = (maximumspace - used)/
                          (2 * sizeof (*fci.tempcodeposforradixsort));
    } else
    {
      fci.buf.allocated /= 2UL;
    }
    if (!onlyaccumulation)
    {
      fci.buf.spaceGtUlongPair = gt_malloc(sizeof (*fci.buf.spaceGtUlongPair)
                                           * fci.buf.allocated);
      fci.tempcodeposforradixsort
        = gt_malloc(sizeof (*fci.tempcodeposforradixsort) * fci.buf.allocated);
    }
    /*
    printf("allocate %.2f for codeposbuffer\n",
             GT_MEGABYTES(sizeof (*fci.buf.spaceGtUlongPair) * 2 *
                          fci.buf.allocated));
    */
    fci.flushcount = 0;
    for (part = 0; !onlyaccumulation &&
                   part < gt_suftabparts_numofparts(suftabparts); part++)
    {
      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "to insert suffixes into buckets",stdout);
      }
      fci.widthofpart = gt_suftabparts_widthofpart(part,suftabparts);
      gt_logger_log(logger,"compute part %u (%.2f%% of all candidates)",part,
                   (double) 100.0 * fci.widthofpart/suftabentries);
      fci.currentminindex = gt_suftabparts_minindex(part,suftabparts);
      fci.currentmaxindex = gt_suftabparts_maxindex(part,suftabparts);
      if (fci.mappedallfirstcodes != NULL)
      {
        fci.tab.allfirstcodes
          = (unsigned long *)
            gt_Sfxmappedrange_map(fci.mappedallfirstcodes,
                                  part,
                                  fci.currentminindex,
                                  fci.currentmaxindex);
      }
      if (fci.mappedcountocc != NULL)
      {
        unsigned long largestindex2map;

        if (part > 0 && fci.tab.overflow_index > 0)
        {
          largestindex2map = MIN(fci.tab.overflow_index-1,fci.currentmaxindex);
        } else
        {
          largestindex2map = fci.currentmaxindex;
        }
        if (fci.currentminindex <= largestindex2map)
        {
          gt_firstcodes_countocc_remap(&fci.tab,
                                       (uint32_t *) gt_Sfxmappedrange_map(
                                                          fci.mappedcountocc,
                                                          part,
                                                          fci.currentminindex,
                                                          largestindex2map));
        }
      }
      if (fci.mappedoverflow != NULL)
      {
        gt_assert(fci.tab.overflow_index > 0);
        if (fci.tab.overflow_index <= fci.currentmaxindex)
        {
          if (fci.currentminindex <= fci.tab.overflow_index)
          {
            fci.tab.overflow_leftborder
              = (unsigned long *)
                gt_Sfxmappedrange_map(fci.mappedoverflow,
                                      part,
                                      0,
                                      fci.currentmaxindex);
          } else
          {
            fci.tab.overflow_leftborder
              = (unsigned long *)
                gt_Sfxmappedrange_map(fci.mappedoverflow,
                                      part,
                                      fci.currentminindex -
                                        fci.tab.overflow_index,
                                      fci.currentmaxindex -
                                        fci.tab.overflow_index);
          }
        }
      }
      if (fci.mappedmarkprefix != NULL)
      {
        gt_marksubstring_bits_map(fci.buf.markprefix,
                                  (GtBitsequence *)
                                  gt_Sfxmappedrange_map(fci.mappedmarkprefix,
                                                        part,
                                                        fci.currentminindex,
                                                        fci.currentmaxindex));
      }
      /*
      gt_logger_log(logger,"current space for part %u=%.2f",
                   part,GT_MEGABYTES(gt_ma_get_space_current() +
                                gt_fa_get_space_current()));*/
      fci.buf.currentmincode = gt_suftabparts_mincode(part,suftabparts);
      fci.buf.currentmaxcode = gt_suftabparts_maxcode(part,suftabparts);
      gt_spmsuftab_partoffset(fci.spmsuftab,
                              gt_suftabparts_offset(part,suftabparts));
      gt_firstcodes_insertsuffix_getencseqkmers_twobitencoding(
                                    encseq,
                                    readmode,
                                    kmersize,
                                    minmatchlength,
                                    &fci.buf,
                                    NULL);
      gt_firstcodes_insertsuffixes_flush(&fci);
      if (part == gt_suftabparts_numofparts(suftabparts) - 1)
      {
        gt_firstcode_delete_before_end(&fci);
      }
      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "to sort buckets of suffixes",stdout);
      }
      if (gt_firstcodes_sortremaining(encseq,
                                  readmode,
                                  srsw,
                                  seqnum_relpos_bucket,
                                  fci.spmsuftab,
                                  fci.buf.snrp,
                                  &fci.tab,
                                  gt_suftabparts_minindex(part,suftabparts),
                                  gt_suftabparts_maxindex(part,suftabparts),
                                  gt_suftabparts_sumofwidth(part,suftabparts),
                                  (unsigned long) kmersize,
                                  itvprocess,
                                  itvprocessdata,
                                  withsuftabcheck,
                                  err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "cleaning up",stdout);
  }
  gt_shortreadsort_delete(srsw);
  gt_free(seqnum_relpos_bucket);
  if (haserr)
  {
    gt_firstcode_delete_before_end(&fci);
    gt_free(fci.buf.spaceGtUlong);
    gt_free(fci.tab.allfirstcodes);
    gt_Sfxmappedrangelist_delete(sfxmrlist);
    gt_free(fci.tempcodeforradixsort);
  }
  gt_suftabparts_delete(suftabparts);
  if (!haserr)
  {
    if (!onlyaccumulation)
    {
      gt_logger_log(logger,"firstcodeposhits=%lu (%.3f%% of all suffixes), "
                           "%u rounds (avg length %lu)",
                           fci.firstcodeposhits,
                           100.0 * (double) fci.firstcodeposhits/totallength,
                           fci.flushcount,
                           fci.firstcodeposhits/fci.flushcount);
      gt_assert(fci.firstcodeposhits == suftabentries);
    } else
    {
      gt_firstcode_delete_before_end(&fci);
    }
  }
  if (haserr || fci.mappedcountocc == NULL)
  {
    gt_firstcodes_countocc_delete(&fci.tab);
  } else
  {
    gt_free(fci.tab.countocc_small);
    gt_free(fci.tab.countocc_samples);
    gt_hashtable_delete(fci.tab.countocc_exceptions);
    gt_str_delete(fci.tab.outfilenameleftborder);
  }
  if (haserr || fci.mappedoverflow == NULL)
  {
    gt_free(fci.tab.overflow_leftborder);
  }
  gt_spmsuftab_delete(fci.spmsuftab);
  gt_Sfxmappedrange_delete(fci.mappedcountocc);
  gt_Sfxmappedrange_delete(fci.mappedoverflow);
  gt_Sfxmappedrange_delete(fci.mappedallfirstcodes);
  gt_seqnumrelpos_delete(fci.buf.snrp);
  if (!haserr)
  {
    gt_log_log("workspace = %.2f",GT_MEGABYTES(fcsl.workspace));
    gt_log_log("size to split = %.2f",GT_MEGABYTES(fcsl.splitspace));
    gt_log_log("suftabsize = %.2f",GT_MEGABYTES(suftab_size));
    gt_log_log("estimatedspace = %.2f",GT_MEGABYTES(fcsl.workspace +
                                                    fcsl.splitspace +
                                                    suftab_size));
  }
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  /*printf("callfirstcodes_find = %lu\n",callfirstcodes_find);*/
  return haserr ? -1 : 0;
}
