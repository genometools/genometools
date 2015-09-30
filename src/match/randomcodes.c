/*
  Copyright (c) 2011-2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011-2013 Center for Bioinformatics, University of Hamburg

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
#include "core/fa.h"
#include "core/arraydef.h"
#include "core/codetype.h"
#include "core/encseq.h"
#include "core/error_api.h"
#include "core/logger_api.h"
#include "core/mathsupport.h"
#include "core/radix_sort.h"
#include "core/showtime.h"
#include "core/spacecalc.h"
#include "core/spacepeak.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/minmax.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif
#include "match/firstcodes-buf.h"
#include "match/firstcodes-spacelog.h"
#include "match/randomcodes-tab.h"
#include "match/firstcodes-accum.h"
#include "match/randomcodes-insert.h"
#include "match/randomcodes.h"
#include "match/seqnumrelpos.h"
#include "match/randomcodes-sfx-partssuf.h"
#include "match/randomcodes-correct.h"
#include "match/sfx-shortreadsort.h"
#include "match/sfx-suffixer.h"
#include "match/spmsuftab.h"
#include "match/kmercodes.h"

typedef struct
{
  GtUword *ptr, code, index;
} GtIndexwithcodeRC;

typedef struct
{
  GtIndexwithcodeRC *spaceGtIndexwithcodeRC;
  GtUword width, nextfreeGtIndexwithcodeRC, allocatedGtIndexwithcodeRC;
  unsigned int depth;
} GtArrayGtIndexwithcodeRC;

typedef struct
{
  GtUword total_count,
                total_inserted,
                countcodes,
                numofcodes,
                codebuffer_total,
                currentminindex,
                currentmaxindex,
                differentcodes, /* a copy of the same value as in tab */
                widthofpart;
  GtArrayGtIndexwithcodeRC binsearchcache;
  unsigned int flushcount,
               bitsforposref;
  GtRadixsortinfo *radixsort_code,
                  *radixsort_codepos;
  GtSpmsuftab *spmsuftab;
  GtSfxmappedrange *mappedleftborder,
                   *mappedallrandomcodes;
  GtUword *allrandomcodes;
  GtFirstcodesspacelog *fcsl;
  GtCodeposbuffer buf;
  GtRandomcodestab tab;
} GtRandomcodesinfo;

static void gt_storerandomcodes(void *processinfo,
                               GT_UNUSED bool firstinrange,
                               GT_UNUSED GtUword pos,
                               GtCodetype code)
{
  GtRandomcodesinfo *fci = (GtRandomcodesinfo *) processinfo;

  gt_assert(fci != NULL && firstinrange &&
            fci->allrandomcodes != NULL &&
            fci->countcodes < fci->numofcodes);
  fci->allrandomcodes[fci->countcodes++] = code;
}

static void gt_randomcodes_fillbinsearchcache(GtRandomcodesinfo *fci,
                                             unsigned int addbscache_depth)
{
  size_t allocbytes = 0;

  fci->binsearchcache.depth
    = addbscache_depth + (unsigned int) log10((double) fci->differentcodes);
  fci->binsearchcache.nextfreeGtIndexwithcodeRC = 0;
  fci->binsearchcache.allocatedGtIndexwithcodeRC
    = 1UL << (fci->binsearchcache.depth+1);
  fci->binsearchcache.width
    = fci->differentcodes/fci->binsearchcache.allocatedGtIndexwithcodeRC;
  if (fci->binsearchcache.allocatedGtIndexwithcodeRC < fci->differentcodes)
  {
    GtUword idx, current = fci->binsearchcache.width;

    allocbytes = sizeof (*fci->binsearchcache.spaceGtIndexwithcodeRC)
                 * fci->binsearchcache.allocatedGtIndexwithcodeRC;
    fci->binsearchcache.spaceGtIndexwithcodeRC = gt_malloc(allocbytes);
    for (idx=0; idx < fci->binsearchcache.allocatedGtIndexwithcodeRC; idx++)
    {
      gt_assert(current < fci->differentcodes);
      fci->binsearchcache.
           spaceGtIndexwithcodeRC[fci->binsearchcache.nextfreeGtIndexwithcodeRC]
                               .ptr = fci->allrandomcodes + current;
      fci->binsearchcache.
           spaceGtIndexwithcodeRC[fci->binsearchcache.nextfreeGtIndexwithcodeRC]
                               .index = current;
      fci->binsearchcache.
        spaceGtIndexwithcodeRC[fci->binsearchcache.nextfreeGtIndexwithcodeRC++]
        .code = fci->allrandomcodes[current];
      current += fci->binsearchcache.width;
    }
  }
  gt_log_log("binsearchcache.depth=%u => " GT_WU " bytes",
             fci->binsearchcache.depth,
             (GtUword) allocbytes);
  GT_FCI_ADDWORKSPACE(fci->fcsl, "binsearchcache", allocbytes);
}

#define SHOWFOUND(F)\
        if ((F) == NULL)\
        {\
          fprintf(stderr, "%s = NULL\n", #F);\
        } else\
        {\
          fprintf(stderr, "%s = " GT_WU "\n", #F, \
                   (GtUword) ((F) - fci->allrandomcodes));\
        }

GtUword gt_randomcodes_find_accu(const GtRandomcodesinfo *fci,
                                       GtUword code)
{
  const GtUword *found = NULL, *leftptr = NULL, *midptr, *rightptr = NULL;

  if (code <= fci->allrandomcodes[0])
  {
    return 0;
  }
  if (fci->binsearchcache.spaceGtIndexwithcodeRC != NULL)
  {
    const GtIndexwithcodeRC *leftic, *midic, *rightic;
    unsigned int depth;

    leftic = fci->binsearchcache.spaceGtIndexwithcodeRC;
    rightic = fci->binsearchcache.spaceGtIndexwithcodeRC +
              fci->binsearchcache.nextfreeGtIndexwithcodeRC - 1;
    for (depth = 0; /* Nothing */; depth++)
    {
      midic = leftic + GT_DIV2((GtUword) (rightic-leftic));
      if (code < midic->code)
      {
        found = midic->ptr;
        if (depth < fci->binsearchcache.depth)
        {
          rightic = midic - 1;
        } else
        {
          gt_assert(leftic->ptr != NULL && rightic->ptr != NULL);
          if (leftic > fci->binsearchcache.spaceGtIndexwithcodeRC)
          {
            leftptr = (leftic-1)->ptr + 1;
          } else
          {
            gt_assert(code > fci->allrandomcodes[0]);
            leftptr = fci->allrandomcodes + 1;
          }
          rightptr = rightic->ptr - 1;
          break;
        }
      } else
      {
        if (code > midic->code)
        {
          if (depth < fci->binsearchcache.depth)
          {
            leftic = midic + 1;
          } else
          {
            gt_assert(leftic->ptr != NULL && rightic->ptr != NULL);
            leftptr = leftic->ptr + 1;
            if (rightic < fci->binsearchcache.spaceGtIndexwithcodeRC +
                          fci->binsearchcache.nextfreeGtIndexwithcodeRC - 1)
            {
              rightptr = (rightic+1)->ptr - 1;
            } else
            {
              rightptr = fci->allrandomcodes + fci->differentcodes - 1;
            }
            break;
          }
        } else
        {
          gt_assert(midic->ptr != NULL);
          return (GtUword) (midic->ptr - fci->allrandomcodes);
        }
      }
    }
    gt_assert(leftptr != NULL && rightptr != NULL);
  } else
  {
    leftptr = fci->allrandomcodes + 1;
    rightptr = fci->allrandomcodes + fci->differentcodes - 1;
  }
  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((GtUword) (rightptr-leftptr));
    if (code < *midptr)
    {
      rightptr = midptr - 1;
      if (code > *rightptr)
      {
        return (GtUword) (midptr - fci->allrandomcodes);
      }
      found = midptr;
    } else
    {
      if (code > *midptr)
      {
        leftptr = midptr + 1;
        if (code <= *leftptr)
        {
          return (GtUword) (leftptr - fci->allrandomcodes);
        }
      } else
      {
        gt_assert(midptr != NULL);
        return (GtUword) (midptr - fci->allrandomcodes);
      }
    }
  }
 return (found != NULL) ? (GtUword) (found - fci->allrandomcodes)
                         : ULONG_MAX;
}

const GtUword *gt_randomcodes_find_insert(const GtRandomcodesinfo *fci,
                                               GtUword code)
{
  const GtUword *found = NULL, *leftptr = NULL, *midptr, *rightptr = NULL;

  leftptr = fci->allrandomcodes + fci->currentminindex;
  rightptr = fci->allrandomcodes + fci->currentmaxindex;
  if (code < *leftptr)
  {
    return leftptr;
  }
  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((GtUword) (rightptr-leftptr));
    if (code < *midptr)
    {
      rightptr = midptr - 1;
      if (code > *rightptr)
      {
        return midptr;
      }
      found = midptr;
    } else
    {
      if (code > *midptr)
      {
        leftptr = midptr + 1;
        if (code <= *leftptr)
        {
          return leftptr;
        }
      } else
      {
        return midptr;
      }
    }
  }
  return found;
}

static GtUword gt_randomcodes_accumulatecounts_merge(
                                        GtRandomcodesinfo *fci,
                                        const GtUword *querystream_fst,
                                        const GtUword *subjectstream_fst)
{
  GtUword found = 0;
  const GtUword *query = querystream_fst,
                      *subject = subjectstream_fst,
                      *querystream_lst = fci->buf.spaceGtUword
                                         + fci->buf.nextfree - 1,
                      *subjectstream_lst = fci->allrandomcodes
                                           + fci->differentcodes - 1;

  while (query <= querystream_lst && subject <= subjectstream_lst)
  {
    if (*query <= *subject)
    {
      gt_randomcodes_countocc_increment(&fci->tab, (GtUword)
          (subject - fci->allrandomcodes));
      found++;
      query++;
    } else
    {
      subject++;
    }
  }
  return found;
}

static void gt_randomcodes_accumulatecounts_flush(void *data)
{
  GtRandomcodesinfo *fci = (GtRandomcodesinfo *) data;

  if (fci->buf.nextfree > 0)
  {
    GtUword foundindex;

    gt_assert(fci->allrandomcodes != NULL);
    fci->codebuffer_total += fci->buf.nextfree;
    gt_radixsort_inplace_sort(fci->radixsort_code, fci->buf.nextfree);
    foundindex = gt_randomcodes_find_accu(fci, fci->buf.spaceGtUword[0]);
    gt_assert(foundindex != ULONG_MAX);
    fci->total_count += gt_randomcodes_accumulatecounts_merge(fci,
        fci->buf.spaceGtUword, fci->allrandomcodes + foundindex);
    gt_assert(fci->total_count == fci->codebuffer_total);
    fci->flushcount++;
    fci->buf.nextfree = 0;
  }
}

static GtUword gt_randomcodes_insertsuffixes_merge(
                                        GtRandomcodesinfo *fci,
                                        const GtUwordPair *querystream_fst,
                                        const GtUword *subjectstream_fst)
{
  GtUword found = 0, idx;
  const GtUwordPair *query = querystream_fst,
                    *querystream_lst = fci->buf.spaceGtUwordPair +
                                       fci->buf.nextfree - 1;
  const GtUword *subject = subjectstream_fst,
                      *subjectstream_lst = fci->allrandomcodes +
                                           fci->currentmaxindex;

  while (query <= querystream_lst && subject <= subjectstream_lst)
  {
    if (query->a <= *subject)
    {
      idx = gt_randomcodes_insertionindex(&fci->tab,
          (GtUword)(subject - fci->allrandomcodes));
      gt_spmsuftab_set(fci->spmsuftab, idx,
          gt_spmsuftab_usebitsforpositions(fci->spmsuftab)
          ? gt_seqnumrelpos_decode_pos(fci->buf.snrp, query->b)
          : query->b);
      found++;
      query++;
    } else
    {
      subject++;
    }
  }
  return found;
}

static void gt_randomcodes_insertsuffixes_flush(void *data)
{
  GtRandomcodesinfo *fci = (GtRandomcodesinfo *) data;

  if (fci->buf.nextfree > 0)
  {
    const GtUword *ptr;

    gt_assert(fci->allrandomcodes != NULL);
    fci->codebuffer_total += fci->buf.nextfree;
    gt_radixsort_inplace_sort(fci->radixsort_codepos, fci->buf.nextfree);
    ptr = gt_randomcodes_find_insert(fci, fci->buf.spaceGtUwordPair[0].a);
    gt_assert(ptr != NULL);
    fci->total_inserted += gt_randomcodes_insertsuffixes_merge(fci,
        fci->buf.spaceGtUwordPair, ptr);
    gt_assert(fci->total_inserted == fci->codebuffer_total);
    fci->flushcount++;
    fci->buf.nextfree = 0;
  }
}

static void gt_randomcodes_checksuftab_bucket(const GtEncseq *encseq,
                                             GtReadmode readmode,
                                             GtEncseqReader *esr1,
                                             GtEncseqReader *esr2,
                                             GtUword previoussuffix,
                                             bool previousdefined,
                                             const GtUword
                                               *seqnum_relpos_bucket,
                                             const GtSeqnumrelpos *snrp,
                                             GT_UNUSED const uint16_t
                                               *lcptab_bucket,
                                             GtUword numberofsuffixes)
{
  GtUword idx, current, maxlcp,
                totallength = gt_encseq_total_length(encseq);
  const GtUword depth = 0;
  GT_UNUSED int cmp;
  const bool specialsareequal = false, specialsareequalatdepth0 = false;

  gt_assert(!previousdefined || previoussuffix < totallength);
  for (idx = 0; idx < numberofsuffixes; idx++)
  {
    current = gt_seqnumrelpos_decode_pos(snrp, seqnum_relpos_bucket[idx]);
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
      gt_assert(idx == 0 || maxlcp == (GtUword) lcptab_bucket[idx]);
    }
    previoussuffix = current;
    previousdefined = true;
  }
}

/* this can be maybe improved by using optimized code for
 * determining the position of the most significant bit */
static GtUword gt_randomcodes_codelcp(GtUword a, GtUword b)
{
  GtUword lcpvalue = 0;
  GtUword mask;
  const GtUword xorvalue = a ^ b;
  for (mask = (GtUword) 3 << GT_MULT2(GT_UNITSIN2BITENC - 1);
      mask > 0; mask >>= 2, lcpvalue++)
    if (xorvalue & mask)
      return lcpvalue;
  return lcpvalue;
}

static int gt_randomcodes_sortremaining(GtShortreadsortworkinfo *srsw,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       const GtSpmsuftab *spmsuftab,
                                       const GtSeqnumrelpos *snrp,
                                       const GtRandomcodestab *rct,
                                       const GtUword *allrandomcodes,
                                       GtUword minindex,
                                       GtUword maxindex,
                                       GtUword sumofwidth,
                                       unsigned int sortingdepth,
                                       unsigned int bucketkeysize,
                                       GtRandomcodesintervalprocess itvprocess,
                                       GtRandomcodesintervalprocess_end
                                              itvprocess_end,
                                       void *itvprocessdata,
                                       bool withsuftabcheck,
                                       GtError *err)
{
  GtUword current,
                next = GT_UNDEF_UWORD,
                idx,
                width,
                sumwidth = 0,
                previoussuffix = 0;
  GtShortreadsortresult srsresult;
  bool previousdefined = false, haserr = false;

  gt_assert(allrandomcodes != NULL);
  current = gt_randomcodes_get_leftborder(rct, minindex);
  for (idx = minindex; idx <= maxindex; idx++)
  {
    if (idx < maxindex)
    {
      next = gt_randomcodes_get_leftborder(rct, idx+1);
      gt_assert(current <= next);
      width = next - current;
    } else
    {
      gt_assert(sumofwidth >= current);
      width = sumofwidth - current;
    }
    sumwidth += width;
    gt_assert(sumwidth <= spmsuftab->numofentries);
    if (width >= 2UL)
    {
      GtUword lcpvalue = 0;
      if (idx > minindex)
      {
        lcpvalue = gt_randomcodes_codelcp(allrandomcodes[idx - 1],
            allrandomcodes[idx]);
        gt_assert(lcpvalue >=
            (GtUword)(GT_UNITSIN2BITENC - bucketkeysize));
        lcpvalue -= (GT_UNITSIN2BITENC - (GtUword) bucketkeysize);
      }
      gt_assert(lcpvalue < (GtUword)GT_UNITSIN2BITENC);
      gt_shortreadsort_firstcodes_sort(&srsresult,
                                       srsw,
                                       snrp,
                                       encseq,
                                       spmsuftab,
                                       current,
                                       width,
                                       lcpvalue,
                                       (GtUword)sortingdepth);
      if (withsuftabcheck)
      {
        gt_randomcodes_checksuftab_bucket(encseq,
                                         readmode,
                                         NULL,
                                         NULL,
                                         previoussuffix,
                                         previousdefined,
                                         srsresult.suftab_bucket,
                                         snrp,
                                         srsresult.lcptab_bucket,
                                         width);
        previousdefined = true;
        previoussuffix
          = gt_seqnumrelpos_decode_pos(snrp, srsresult.suftab_bucket[width-1]);
      }
      if (itvprocess != NULL)
      {
        if (itvprocess(itvprocessdata,
                       srsresult.suftab_bucket,
                       snrp, srsresult.lcptab_bucket, width,
                       sortingdepth, err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
    else if (width == 1UL)
    {
      if (itvprocess != NULL)
      {
        GtUword suftab_value = gt_spmsuftab_get(spmsuftab, current);
        const uint16_t lcptab_value = 0;
        if (gt_spmsuftab_usebitsforpositions(spmsuftab))
        {
          GtUword seqnum, relpos;
          seqnum = gt_encseq_seqnum(encseq, suftab_value);
          relpos = suftab_value - gt_encseq_seqstartpos(encseq, seqnum);
          suftab_value = gt_seqnumrelpos_encode(snrp, seqnum, relpos);
        }
        if (itvprocess(itvprocessdata, &suftab_value, snrp, &lcptab_value, 1UL,
                       sortingdepth, err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
    gt_assert((minindex == maxindex) || next != GT_UNDEF_UWORD);
    current = next;
  }
  if (itvprocess_end != NULL)
  {
    itvprocess_end(itvprocessdata);
  }
  return haserr ? -1 : 0;
}

#ifdef GT_THREADS_ENABLED

static GtUword gt_randomcodes_findfirstlarger(const GtRandomcodestab *rct,
                                                   GtUword start,
                                                   GtUword end,
                                                   GtUword offset)
{
  GtUword left = start, right = end, found = end, mid, midval;

  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_randomcodes_get_leftborder(rct, mid);
    if (offset == midval)
    {
      return mid;
    }
    if (offset < midval)
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

static GtUword *gt_randomcodes_evenly_divide_part(
    const GtRandomcodestab *rct, GtUword partminindex,
    GtUword partmaxindex, GtUword numofsuffixes,
    unsigned int numofparts)
{
  GtUword *endindexes, widthofpart, offset;
  unsigned int part, remainder;

  gt_assert(partminindex < partmaxindex && numofparts >= 2U);
  widthofpart = numofsuffixes/numofparts;
  endindexes = gt_malloc(sizeof (*endindexes) * numofparts);
  offset = gt_randomcodes_get_leftborder(rct, partminindex);
  remainder = (unsigned int) (numofsuffixes % (GtUword) numofparts);
  for (part=0; part < numofparts; part++)
  {
    if (remainder > 0)
    {
      offset += widthofpart + 1;
      remainder--;
    } else
    {
      offset += widthofpart;
    }
    if (part == numofparts - 1)
    {
      endindexes[part] = partmaxindex;
    } else
    {
      GtUword start = (part == 0) ? partminindex : endindexes[part-1] + 1;

      endindexes[part] = gt_randomcodes_findfirstlarger(rct, start,
                                                        partmaxindex, offset);
      gt_assert(endindexes[part] <= partmaxindex);
    }
  }
  return endindexes;
}

typedef struct
{
  GtShortreadsortworkinfo *srsw;
  const GtEncseq *encseq;
  GtReadmode readmode;
  const GtSpmsuftab *spmsuftab;
  const GtSeqnumrelpos *snrp;
  const GtRandomcodestab *rct;
  const GtUword *allrandomcodes;
  GtUword minindex,
                maxindex,
                sumofwidth;
  unsigned int sortingdepth,
               bucketkeysize;
  bool withsuftabcheck;
  GtRandomcodesintervalprocess itvprocess;
  GtRandomcodesintervalprocess_end itvprocess_end;
  void *itvprocessdata;
  GtError *err;
  GtThread *thread;
} GtRandomcodesSortRemainingThreadinfo;

static void *gt_randomcodes_thread_caller_sortremaining(void *data)
{
  GtRandomcodesSortRemainingThreadinfo *threadinfo =
    (GtRandomcodesSortRemainingThreadinfo *) data;

  if (gt_randomcodes_sortremaining(threadinfo->srsw,
                                  threadinfo->encseq,
                                  threadinfo->readmode,
                                  threadinfo->spmsuftab,
                                  threadinfo->snrp,
                                  threadinfo->rct,
                                  threadinfo->allrandomcodes,
                                  threadinfo->minindex,
                                  threadinfo->maxindex,
                                  threadinfo->sumofwidth,
                                  threadinfo->sortingdepth,
                                  threadinfo->bucketkeysize,
                                  threadinfo->itvprocess,
                                  threadinfo->itvprocess_end,
                                  threadinfo->itvprocessdata,
                                  threadinfo->withsuftabcheck,
                                  threadinfo->err) != 0)
  {
    gt_assert(false);
  }
  return NULL;
}

static int gt_randomcodes_thread_sortremaining(
                                       GtShortreadsortworkinfo **srswtab,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       const GtSpmsuftab *spmsuftab,
                                       const GtSeqnumrelpos *snrp,
                                       const GtRandomcodestab *rct,
                                       GtUword *allrandomcodes,
                                       GtUword partminindex,
                                       GtUword partmaxindex,
                                       GtUword widthofpart,
                                       GtUword sumofwidth,
                                       unsigned int sortingdepth,
                                       unsigned int bucketkeysize,
                                       GtRandomcodesintervalprocess itvprocess,
                                       GtRandomcodesintervalprocess_end
                                         itvprocess_end,
                                       void *itvprocessdatatab,
                                       bool withsuftabcheck,
                                       unsigned int threads,
                                       GtLogger *logger,
                                       GtError *err)
{
  unsigned int t;
  GtUword sum = 0, *endindexes;
  GtRandomcodesSortRemainingThreadinfo *threadinfo;
  bool haserr = false;

  gt_assert(threads >= 2U);
  endindexes = gt_randomcodes_evenly_divide_part(rct, partminindex,
      partmaxindex, widthofpart, threads);
  threadinfo = gt_malloc(sizeof (*threadinfo) * threads);
  for (t=0; t<threads; t++)
  {
    GtUword lb;

    threadinfo[t].srsw = srswtab[t];
    threadinfo[t].encseq = encseq;
    threadinfo[t].spmsuftab = spmsuftab;
    threadinfo[t].snrp = snrp;
    threadinfo[t].rct = rct;
    threadinfo[t].allrandomcodes = allrandomcodes;
    threadinfo[t].minindex = t == 0 ? partminindex : endindexes[t-1] + 1,
    threadinfo[t].maxindex = endindexes[t];
    threadinfo[t].readmode = readmode;
    threadinfo[t].withsuftabcheck = withsuftabcheck;
    threadinfo[t].sortingdepth = sortingdepth;
    threadinfo[t].bucketkeysize = bucketkeysize;
    threadinfo[t].itvprocess = itvprocess;
    threadinfo[t].itvprocess_end = itvprocess_end;
    if (itvprocessdatatab == NULL)
    {
      threadinfo[t].itvprocessdata = NULL;
    } else
    {
      threadinfo[t].itvprocessdata = ((void **) itvprocessdatatab)[t];
    }
    threadinfo[t].err = err;
    lb = gt_randomcodes_get_leftborder(rct, threadinfo[t].minindex);
    if (t < threads - 1)
    {
      threadinfo[t].sumofwidth
        = gt_randomcodes_get_leftborder(rct, threadinfo[t].maxindex+1);
    } else
    {
      threadinfo[t].sumofwidth = sumofwidth;
    }
    gt_assert(lb < threadinfo[t].sumofwidth);
    gt_logger_log(logger, "thread %u: process [" GT_WU ", " GT_WU "]=[" GT_WU
                  ", " GT_WU "] of width " GT_WU,
                  t, threadinfo[t].minindex, threadinfo[t].maxindex, lb,
                  threadinfo[t].sumofwidth, threadinfo[t].sumofwidth - lb);
    sum += threadinfo[t].sumofwidth - lb;
    threadinfo[t].thread
      = gt_thread_new (gt_randomcodes_thread_caller_sortremaining,
                       threadinfo + t, err);
    if (threadinfo[t].thread == NULL)
    {
      haserr = true;
    }
  }
  gt_assert (haserr || sum == widthofpart);
  for (t=0; t<threads; t++)
  {
    if (!haserr)
    {
      gt_thread_join(threadinfo[t].thread);
    }
    gt_thread_delete(threadinfo[t].thread);
  }
  gt_free(threadinfo);
  gt_free(endindexes);
  return haserr ? -1 : 0;
}
#endif

static void gt_randomcodes_delete_before_end(GtRandomcodesinfo *fci)
{
  if (fci->binsearchcache.spaceGtIndexwithcodeRC != NULL)
  {
    GT_FCI_SUBTRACTWORKSPACE(fci->fcsl, "binsearchcache");
    GT_FREEARRAY(&fci->binsearchcache, GtIndexwithcodeRC);
  }
  if (fci->radixsort_codepos != NULL)
  {
    gt_radixsort_delete(fci->radixsort_codepos);
    GT_FCI_SUBTRACTWORKSPACE(fci->fcsl, "radixsort_codepos");
    fci->radixsort_codepos = NULL;
  }
}

static GtUword gt_randomcodes_idx2mincode(const GtRandomcodesinfo *fci,
                                            GtUword idx)
{
  gt_assert(idx <= fci->differentcodes);
  if (idx == 0)
  {
    return 0;
  }
  else
  {
    return fci->allrandomcodes[idx-1];
  }
}

static GtUword gt_randomcodes_idx2maxcode(const GtRandomcodesinfo *fci,
                                            GtUword idx)
{
  gt_assert(idx <= fci->differentcodes);
  if (idx == fci->differentcodes)
  {
    return fci->allrandomcodes[idx-1];
  }
  return fci->allrandomcodes[idx];
}

void gt_rungetencseqkmers_rc(const GtEncseq *encseq, unsigned int bucketkeysize)
{
  const GtReadmode readmode = GT_READMODE_FORWARD;

  getencseqkmers_twobitencoding(encseq,
                                readmode,
                                bucketkeysize,
                                bucketkeysize,
                                false,
                                NULL,
                                NULL,
                                NULL,
                                NULL);
}

static void run_allcodes_distribution(const GtUword *allrandomcodes,
                                      GtUword differentcodes)
{
  GtUword idx, diff, mindiff = 0, maxdiff = 0, distbits[64+1] = {0};

  for (idx = 1UL; idx < differentcodes; idx++)
  {
    gt_assert(allrandomcodes[idx-1] < allrandomcodes[idx]);
    diff = allrandomcodes[idx] - allrandomcodes[idx-1];
    if (idx == 1UL || diff < mindiff)
    {
      mindiff = diff;
    }
    if (diff > maxdiff)
    {
      maxdiff = diff;
    }
    distbits[gt_determinebitspervalue(diff)]++;
  }
  printf("allrandomcodes: mindiff=" GT_WU ", maxdiff=" GT_WU "(%u bits)\n",
         mindiff, maxdiff, gt_determinebitspervalue(maxdiff));
  for (idx = 0; idx <= 64UL; idx++)
  {
    if (distbits[idx] > 0)
    {
      printf("" GT_WU " bits: " GT_WU "\n", idx, distbits[idx]);
    }
  }
}

static int gt_randomcodes_init(GtRandomcodesinfo *fci,
                              const GtEncseq *encseq,
                              bool withsuftabcheck,
                              unsigned int skipshorter,
                              GtError *err)
{
  GtUword maxseqlength, maxrelpos, numofsequences;
  unsigned int bitsforrelpos, bitsforseqnum;
  bool haserr = false;

  maxseqlength = gt_encseq_max_seq_length(encseq);
  if (maxseqlength > (GtUword) skipshorter)
  {
    maxrelpos = maxseqlength - (GtUword) skipshorter;
  } else
  {
    maxrelpos = 0;
  }
  bitsforrelpos = gt_determinebitspervalue(maxrelpos);
  fci->buf.snrp = gt_seqnumrelpos_new(bitsforrelpos, encseq);
  fci->buf.markprefix = NULL;
  fci->buf.marksuffix = NULL;
  fci->buf.accum_all = true;
  numofsequences = gt_encseq_num_of_sequences(encseq);
  gt_assert(numofsequences > 0);
  bitsforseqnum = gt_determinebitspervalue(numofsequences - 1);
  if (bitsforseqnum + bitsforrelpos > (unsigned int) GT_INTWORDSIZE)
  {
    gt_seqnumrelpos_delete(fci->buf.snrp);
    fci->buf.snrp = NULL;
    gt_error_set(err, "cannot process encoded sequences with " GT_WU
                 " sequences of length up to " GT_WU " (%u+%u bits)",
                 numofsequences, maxseqlength, bitsforseqnum, bitsforrelpos);
    haserr = true;
  }
  fci->fcsl = gt_firstcodes_spacelog_new();
  fci->spmsuftab = NULL;
  fci->radixsort_code = NULL;
  fci->radixsort_codepos = NULL;
  fci->buf.spaceGtUwordPair = NULL;
  fci->buf.spaceGtUword = NULL;
  fci->mappedallrandomcodes = NULL;
  fci->mappedleftborder = NULL;
  GT_FCI_ADDWORKSPACE(fci->fcsl, "encseq", (size_t)
                      gt_encseq_sizeofrep(encseq));
  if (withsuftabcheck)
  {
    gt_firstcodes_spacelog_start_diff(fci->fcsl);
  }
  fci->allrandomcodes = NULL;
  gt_randomcodes_countocc_setnull(&fci->tab);
  fci->countcodes = 0;
  fci->total_count = 0;
  fci->total_inserted = 0;
  fci->bitsforposref = bitsforseqnum + bitsforrelpos;
  GT_INITARRAY(&fci->binsearchcache, GtIndexwithcodeRC);
  return haserr ? -1 : 0;
}

#define GT_RANDOMCODES_NOFSAMPLES_MIN 2UL

static inline GtUword gt_randomcodes_calculate_nofsamples(
    GT_UNUSED const GtEncseq *encseq, GtUword nofsequences,
    GtUword totallength, unsigned int bucketkeysize,
    unsigned int sampling_factor)
{
  GtUword nofkmers = totallength,
                nofnonkmers = (bucketkeysize + 1) * nofsequences - 1,
                nofsamples;
  gt_log_log("totallength = " GT_WU "", nofkmers);
  gt_log_log("nofsequences = " GT_WU "", nofsequences);
  if (nofnonkmers < nofkmers)
  {
    nofkmers -= nofnonkmers;
    gt_log_log("nofkmers = " GT_WU "", nofkmers);
  }
  else
  {
    gt_assert(gt_encseq_min_seq_length(encseq) <= (GtUword)bucketkeysize);
  }
  nofsamples = nofkmers / sampling_factor;
  return MAX(GT_RANDOMCODES_NOFSAMPLES_MIN, nofsamples);
}

static void gt_randomcodes_generate_sampling_positions(GtUword *buffer,
    GtUword numofsamples, const GtEncseq *encseq,
    GtUword totallength, unsigned int sampling_factor,
    unsigned int bucketkeysize, bool sort, GtTimer *timer)
{
  GtUword i, randmax;
  bool sorted = true;
  randmax = (GtUword)(GT_MULT2(sampling_factor) -
      GT_DIV16(sampling_factor));
  if (randmax >= totallength)
    randmax = totallength;
  buffer[0] = gt_rand_max(randmax);
  for (i = 1UL; i < numofsamples; i++)
  {
    GtUword sn, sp, sl, rp;
    buffer[i] = buffer[i-1];
    while (true)
    {
      buffer[i] += (gt_rand_max(randmax) + 1UL);
      if (buffer[i] >= totallength)
      {
        buffer[i] = 0;
        sorted = false;
      }
      if (!gt_encseq_position_is_separator(encseq, buffer[i],
            GT_READMODE_FORWARD))
      {
        sn = gt_encseq_seqnum(encseq, buffer[i]);
        sp = gt_encseq_seqstartpos(encseq, sn);
        sl = gt_encseq_seqlength(encseq, sn);
        rp = buffer[i] - sp;
        if (rp < sl - bucketkeysize)
          break;
      }
    }
  }
  if (sort && !sorted)
  {
    if (timer != NULL)
      gt_timer_show_progress(timer, "to sort sampling positions", stdout);
    gt_radixsort_inplace_ulong(buffer, numofsamples);
  }
  if (sorted)
  {
    gt_log_log("range of sampling positions = [" GT_WU ", " GT_WU "]",
        buffer[0], buffer[numofsamples - 1]);
  }
}

static int gt_randomcodes_collectcodes(GtRandomcodesinfo *fci,
    bool usefirstcodes, unsigned int sampling_factor, const GtEncseq *encseq,
    GtReadmode readmode, unsigned int bucketkeysize, size_t maximumspace,
    GtLogger *logger, GtTimer *timer, GtError *err)
{
  size_t sizeforcodestable;
  GtUword totallength = gt_encseq_total_length(encseq),
                nofsequences = gt_encseq_num_of_sequences(encseq);
  GtUword maskright = GT_MASKRIGHT((GtUword) bucketkeysize);

  if (usefirstcodes)
  {
    fci->numofcodes = nofsequences + 1;
  }
  else
  {
    fci->numofcodes = gt_randomcodes_calculate_nofsamples(encseq,
        nofsequences, totallength, bucketkeysize, sampling_factor) + 1;
  }
  sizeforcodestable = sizeof (*fci->allrandomcodes) * fci->numofcodes;
  if (maximumspace > 0 &&
      gt_firstcodes_spacelog_total(fci->fcsl) + sizeforcodestable >
      maximumspace)
  {
    gt_error_set(err, "already used %.2f MB of memory and require %.2f for the "
        "codes table => cannot compute index in at most %.2f MB",
        GT_MEGABYTES(gt_firstcodes_spacelog_total(fci->fcsl)),
        GT_MEGABYTES(sizeforcodestable),
        GT_MEGABYTES(maximumspace));
    return -1;
  }
  fci->allrandomcodes = gt_malloc(sizeforcodestable);
  GT_FCI_ADDSPLITSPACE(fci->fcsl, "allrandomcodes", sizeforcodestable);

  if (usefirstcodes)
  {
    if (timer != NULL)
      gt_timer_show_progress(timer, "to collect first codes", stdout);
    getencseqkmers_twobitencoding(encseq, readmode, bucketkeysize,
        bucketkeysize, true, gt_storerandomcodes, fci, NULL, NULL);
  }
  else
  {
    GtUword i;
    const GtTwobitencoding *twobitenc =
        gt_encseq_twobitencoding_export(encseq);
    GtUword realtotallength = totallength;
    if (gt_encseq_is_mirrored(encseq))
    {
      realtotallength >>= 1;
    }
    if (timer != NULL)
      gt_timer_show_progress(timer, "to generate sampling positions", stdout);
    gt_randomcodes_generate_sampling_positions(fci->allrandomcodes,
        fci->numofcodes - 1, encseq, totallength, sampling_factor,
        bucketkeysize, true, timer);
    if (timer != NULL)
      gt_timer_show_progress(timer, "to collect sample codes", stdout);
    for (i = 0; i < fci->numofcodes - 1; i++)
    {
      GtUword pos = fci->allrandomcodes[i];
      bool revcompl = false;
      if (pos > realtotallength)
      {
        pos = GT_REVERSEPOS(totallength, pos);
        revcompl = true;
      }
      fci->allrandomcodes[i] = gt_kmercode_at_position(twobitenc,
          pos, bucketkeysize);
      if (revcompl)
      {
        fci->allrandomcodes[i] = gt_kmercode_complement(
         gt_kmercode_reverse(fci->allrandomcodes[i], bucketkeysize), maskright);
      }
    }
    fci->countcodes = fci->numofcodes - 1;
  }
  /* add an artificial last bucket to collect suffixes
   * larger than the last code */
  gt_storerandomcodes(fci, true, /* unused*/ 0, maskright);
  gt_logger_log(logger, "have stored " GT_WU " bucket keys", fci->countcodes);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "to sort bucket keys", stdout);
  }
  gt_radixsort_inplace_ulong(fci->allrandomcodes, fci->numofcodes);
  gt_assert(fci->allrandomcodes != NULL);
  fci->differentcodes = gt_randomcodes_remdups(fci->allrandomcodes,
      bucketkeysize, fci->numofcodes, logger);
  if (fci->differentcodes > 0 && fci->differentcodes < fci->numofcodes)
  {
    fci->allrandomcodes = gt_realloc(fci->allrandomcodes,
                                    sizeof (*fci->allrandomcodes) *
                                    fci->differentcodes);
    GT_FCI_SUBTRACTADDSPLITSPACE(fci->fcsl, "allrandomcodes",
                                 sizeof (*fci->allrandomcodes) *
                                 fci->differentcodes);
  }
  gt_randomcodes_countocc_new(fci->fcsl, &fci->tab, fci->differentcodes);
  return 0;
}

static int gt_randomcodes_allocspace(GtRandomcodesinfo *fci,
                                    unsigned int numofparts,
                                    GtUword maximumspace,
                                    GtUword phase2extra,
                                    GtError *err)
{
  if (maximumspace > 0)
  {
    if ((GtUword) gt_firstcodes_spacelog_total(fci->fcsl) +
                        phase2extra >= maximumspace)
    {
      gt_error_set(err, "already used %.2f MB of memory and need %.2f MB later "
                       "=> cannot compute index in at most %.2f MB",
                       GT_MEGABYTES(gt_firstcodes_spacelog_total(fci->fcsl)),
                       GT_MEGABYTES(phase2extra),
                       GT_MEGABYTES(maximumspace));
      return -1;
    } else
    {
      size_t remainspace = (size_t) maximumspace -
                           (gt_firstcodes_spacelog_total(fci->fcsl) +
                            phase2extra);

      fci->buf.allocated
        = gt_radixsort_max_num_of_entries_ulong(remainspace);
      if (fci->buf.allocated < fci->differentcodes / 16UL)
      {
        fci->buf.allocated = fci->differentcodes / 16UL;
      }
    }
  } else
  {
    if (numofparts == 0)
    {
      fci->buf.allocated = gt_radixsort_max_num_of_entries_ulong(
                                 gt_firstcodes_spacelog_total(fci->fcsl)/7UL);
    } else
    {
      fci->buf.allocated = fci->differentcodes/5;
    }
  }
  if (fci->buf.allocated < 16UL)
  {
    fci->buf.allocated = 16UL;
  }
  return 0;
}

static void gt_randomcodes_accumulatecounts_run(GtRandomcodesinfo *fci,
                                               const GtEncseq *encseq,
                                               unsigned int bucketkeysize,
                                               unsigned int skipshorter,
                                               GtLogger *logger,
                                               GtTimer *timer)
{
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "to accumulate counts", stdout);
  }
  gt_assert(fci->buf.allocated > 0);
  fci->radixsort_code = gt_radixsort_new_ulong(fci->buf.allocated);
  fci->buf.spaceGtUword = gt_radixsort_space_ulong(fci->radixsort_code);
  GT_FCI_ADDWORKSPACE(fci->fcsl, "radixsort_code",
                      gt_radixsort_size(fci->radixsort_code));
  fci->buf.fciptr = fci; /* as we need to give fci to the flush function */
  fci->buf.flush_function = gt_randomcodes_accumulatecounts_flush;
  gt_logger_log(logger, "maximum space for accumulating counts %.2f MB",
                GT_MEGABYTES(gt_firstcodes_spacelog_total(fci->fcsl)));
  gt_firstcodes_accum_runkmerscan(encseq, bucketkeysize, skipshorter,
      &fci->buf);
  gt_randomcodes_accumulatecounts_flush(fci);
  gt_logger_log(logger, "codebuffer_total=" GT_WU " (%.3f%% of all suffixes)",
                fci->codebuffer_total,
                100.0 * (double) fci->codebuffer_total/
                                 gt_encseq_total_length(encseq));
  gt_assert(fci->total_count == fci->codebuffer_total);
  if (fci->total_count > 0)
  {
    gt_assert(fci->flushcount > 0);
    gt_logger_log(logger, "total count=" GT_WU " (%.3f%% of all suffixes), "
                         "%u rounds (avg length " GT_WU ")",
                         fci->total_count,
                         100.0 * (double) fci->total_count/
                                          gt_encseq_total_length(encseq),
                         fci->flushcount,
                         fci->codebuffer_total/fci->flushcount);
  }
  gt_radixsort_delete(fci->radixsort_code);
  fci->radixsort_code = NULL;
  GT_FCI_SUBTRACTWORKSPACE(fci->fcsl, "radixsort_code");
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "to compute partial sums", stdout);
  }
}

static void gt_randomcodes_map_sections(GtRandomcodesinfo *fci,
                                       GtSfxmappedrangelist *sfxmrlist)
{
  if (fci->differentcodes > 0)
  {
    fci->mappedallrandomcodes = gt_Sfxmappedrange_new("allrandomcodes",
                                                    fci->differentcodes,
                                                    GtSfxunsignedlong,
                                                    NULL,
                                                    NULL);
    gt_Sfxmappedrangelist_add(sfxmrlist, fci->mappedallrandomcodes);
    fci->mappedleftborder = gt_Sfxmappedrange_new("leftborder",
                                                 fci->differentcodes+1,
                                                 GtSfxuint32_t,
                                                 NULL,
                                                 NULL);
    gt_Sfxmappedrangelist_add(sfxmrlist, fci->mappedleftborder);
  }
}

static int gt_randomcodes_auto_parts(GtRandomcodesinfo *fci,
                                    GtSfxmappedrangelist *sfxmrlist,
                                    unsigned int numofparts,
                                    GtUword *maximumspace,
                                    GtUword maxbucketsize,
                                    unsigned int sortingdepth,
                                    GtUword totallength,
                                    GtUword maxseqlength,
                                    GtUword suftabentries,
                                    GtUword phase2extra,
                                    GtError *err)
{
  int retval;
  GtUword leftbordersize_all;

  if (numofparts == 0 && *maximumspace == 0)
  {
    *maximumspace = (GtUword)
                    (gt_firstcodes_spacelog_peak(fci->fcsl) +
                     phase2extra +
                     gt_shortreadsort_size(true, maxbucketsize,
                       sortingdepth == 0 ? maxseqlength :
                       (GtUword)sortingdepth) + 4UL * 4096UL);
  } else
  {
    gt_assert(*maximumspace > 0);
  }
  if (fci->mappedleftborder != NULL)
  {
    leftbordersize_all = gt_Sfxmappedrange_size_mapped(
                             fci->mappedleftborder,0,
                             gt_randomcodes_leftborder_entries(&fci->tab)-1);
  } else
  {
    leftbordersize_all = 0;
  }
  retval = gt_suftabparts_rc_fit_memlimit(
      gt_firstcodes_spacelog_total(fci->fcsl)
      /*as this is subtracted*/
      + leftbordersize_all
      + phase2extra,
      *maximumspace,
      NULL,
      &fci->tab,
      sfxmrlist,
      totallength,
      fci->bitsforposref,
      0, /* special characters not used */
      suftabentries,
      false, /* suftabuint not used */
      err);
  if (retval < 0)
  {
    return -1;
  } else
  {
    gt_assert(retval > 0);
    return retval;
  }
}

static void gt_randomcodes_allocsize_for_insertion(GtRandomcodesinfo *fci,
                                                  GtUword maximumspace,
                                                  const GtSuftabparts_rc
                                                    *suftabparts_rc,
                                                  GtUword phase2extra)
{
  if (maximumspace > 0)
  {
    const GtUword maxrounds = 400UL;
    size_t used = gt_firstcodes_spacelog_workspace(fci->fcsl) +
                  phase2extra +
                  gt_suftabparts_rc_largestsizemappedpartwise(suftabparts_rc);

    if ((GtUword) used < maximumspace)
    {
      fci->buf.allocated = gt_radixsort_max_num_of_entries_ulongpair(
                                                (size_t) maximumspace - used);
    } else
    {
      fci->buf.allocated /= 4UL;
    }
    if ((GtUword) (fci->codebuffer_total+fci->numofcodes)/
                        fci->buf.allocated > maxrounds)
    {
      fci->buf.allocated
        = (fci->codebuffer_total+fci->numofcodes)/maxrounds;
    }
  } else
  {
    fci->buf.allocated /= 2UL;
  }
}

static int gt_randomcodes_process_part(GtRandomcodesinfo *fci,
                                      const GtEncseq *encseq,
                                      GtReadmode readmode,
                                      unsigned int bucketkeysize,
                                      unsigned int sortingdepth,
                                      unsigned int skipshorter,
                                      const GtSuftabparts_rc *suftabparts_rc,
                                      unsigned int part,
                                      GtUword maximumspace,
#ifndef GT_THREADS_ENABLED
                                      GT_UNUSED
#endif
                                      unsigned int threads,
                                      GtUword suftabentries,
                                      bool withsuftabcheck,
                                      GtShortreadsortworkinfo **srswtab,
                                      GtRandomcodesintervalprocess itvprocess,
                                      GtRandomcodesintervalprocess_end
                                        itvprocess_end,
                                      void *itvprocessdatatab,
                                      GtLogger *logger,
                                      GtTimer *timer,
                                      GtError *err)
{
  GtUword spaceforbucketprocessing = 0;
  void *mapptr;
  bool haserr = false;

  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "to insert suffixes into buckets", stdout);
  }
  fci->widthofpart = gt_suftabparts_rc_widthofpart(part, suftabparts_rc);
  gt_logger_log(logger, "compute part %u (%.2f%% of all kmers)", part,
               (double) 100.0 * fci->widthofpart/suftabentries);
  fci->currentminindex = gt_suftabparts_rc_minindex(part, suftabparts_rc);
  fci->currentmaxindex = gt_suftabparts_rc_maxindex(part, suftabparts_rc);
  if (fci->mappedallrandomcodes != NULL)
  {
    fci->allrandomcodes
      = (GtUword *)
        gt_Sfxmappedrange_map(fci->mappedallrandomcodes,
                              fci->currentminindex,
                              fci->currentmaxindex);
    GT_FCI_ADDSPLITSPACE(fci->fcsl, "allrandomcodes",
                         (size_t) gt_Sfxmappedrange_size_mapped(
                                        fci->mappedallrandomcodes,
                                        fci->currentminindex,
                                        fci->currentmaxindex));
  }
  gt_assert(fci->allrandomcodes != NULL);
  gt_assert(fci->mappedleftborder != NULL);
  mapptr = gt_Sfxmappedrange_map(fci->mappedleftborder,
                                 fci->currentminindex,
                                 fci->currentmaxindex);
  gt_randomcodes_leftborder_remap(&fci->tab, (uint32_t *) mapptr);
  GT_FCI_ADDSPLITSPACE(fci->fcsl, "leftborder",
                       (size_t) gt_Sfxmappedrange_size_mapped(
                                               fci->mappedleftborder,
                                               fci->currentminindex,
                                               fci->currentmaxindex));
  gt_logger_log(logger, "maximum space for part %u: %.2f MB",
                part, GT_MEGABYTES(gt_firstcodes_spacelog_total(fci->fcsl)));
  fci->buf.currentmincode = gt_randomcodes_idx2mincode(fci,
      fci->currentminindex);
  fci->buf.currentmaxcode = gt_randomcodes_idx2maxcode(fci,
      fci->currentmaxindex);
  gt_spmsuftab_partoffset(fci->spmsuftab,
                          gt_suftabparts_rc_offset(part, suftabparts_rc));
  gt_randomcodes_insert_runkmerscan(encseq,
                                   bucketkeysize,
                                   skipshorter,
                                   &fci->buf);
  gt_randomcodes_insertsuffixes_flush(fci);
  if (part == gt_suftabparts_rc_numofparts(suftabparts_rc) - 1)
  {
    gt_randomcodes_delete_before_end(fci);
  }
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "to sort buckets of suffixes", stdout);
  }
  if (maximumspace > 0)
  {
    if ((GtUword) gt_firstcodes_spacelog_total(fci->fcsl)
        < maximumspace)
    {
      spaceforbucketprocessing = maximumspace -
                                 (GtUword)
                                 gt_firstcodes_spacelog_total(fci->fcsl);
      gt_log_log("space left for sortremaining: %.2f",
                 GT_MEGABYTES(spaceforbucketprocessing));
    }
  }
#ifdef GT_THREADS_ENABLED
  if (threads > 1U)
  {
    if (gt_randomcodes_thread_sortremaining(
          srswtab,
          encseq,
          readmode,
          fci->spmsuftab,
          fci->buf.snrp,
          &fci->tab,
          fci->allrandomcodes,
          fci->currentminindex,
          fci->currentmaxindex,
          gt_suftabparts_rc_widthofpart(part, suftabparts_rc),
          gt_suftabparts_rc_sumofwidth(part, suftabparts_rc),
          sortingdepth,
          bucketkeysize,
          itvprocess,
          itvprocess_end,
          itvprocessdatatab,
          withsuftabcheck,
          threads,
          logger,
          NULL) != 0)
    {
      haserr = true;
    }
  } else
#endif
  {
    if (gt_randomcodes_sortremaining(srswtab[0],
          encseq,
          readmode,
          fci->spmsuftab,
          fci->buf.snrp,
          &fci->tab,
          fci->allrandomcodes,
          fci->currentminindex,
          fci->currentmaxindex,
          gt_suftabparts_rc_sumofwidth(part, suftabparts_rc),
          sortingdepth,
          bucketkeysize,
          itvprocess,
          itvprocess_end,
          itvprocessdatatab == NULL
          ? NULL
          : ((void **) itvprocessdatatab)[0],
          withsuftabcheck,
          err) != 0)
    {
      haserr = true;
    }
  }
  if (fci->mappedallrandomcodes != NULL)
  {
    gt_Sfxmappedrange_unmap(fci->mappedallrandomcodes);
    GT_FCI_SUBTRACTSPLITSPACE(fci->fcsl, "allrandomcodes");
  }
  if (fci->mappedleftborder != NULL)
  {
    gt_Sfxmappedrange_unmap(fci->mappedleftborder);
    GT_FCI_SUBTRACTSPLITSPACE(fci->fcsl, "leftborder");
  }
  return haserr ? -1 : 0;
}

int storerandomcodes_getencseqkmers_twobitencoding(const GtEncseq *encseq,
    unsigned int bucketkeysize,
    unsigned int numofparts,
    GtUword maximumspace,
    unsigned int sortingdepth,
    unsigned int skipshorter,
    bool usefirstcodes,
    unsigned int sampling_factor,
    bool withsuftabcheck,
    bool onlyaccumulation,
    bool onlyallrandomcodes,
    GT_UNUSED
    unsigned int addbscache_depth,
    GtUword phase2extra,
    GT_UNUSED bool radixsmall,      /* set to true */
    GT_UNUSED unsigned int radixparts, /* set to 2U */
    GtRandomcodesintervalprocess
    itvprocess,
    GtRandomcodesintervalprocess_end
    itvprocess_end,
    void *itvprocessdatatab,
    GtLogger *logger,
    GtTimer *timer,
    GtError *err)
{
  GtRandomcodesinfo fci;
  size_t suftab_size = 0;
  unsigned int part, threadcount;
  GtUword maxbucketsize, suftabentries = 0, largest_width,
                totallength = gt_encseq_total_length(encseq),
                maxseqlength = gt_encseq_max_seq_length(encseq);
  GtSfxmappedrangelist *sfxmrlist = NULL;
  GtSuftabparts_rc *suftabparts_rc = NULL;
  GtShortreadsortworkinfo **srswtab = NULL;
  const GtReadmode readmode = GT_READMODE_FORWARD;
  bool haserr = false;
#ifdef GT_THREADS_ENABLED
  const unsigned int threads = gt_jobs;
#else
  const unsigned int threads = 1U;
#endif

  if (maxseqlength < (GtUword)skipshorter)
  {
    return 0;
  }
  if (gt_randomcodes_init(&fci,
                         encseq,
                         withsuftabcheck,
                         skipshorter,
                         err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    sfxmrlist = gt_Sfxmappedrangelist_new();
    if (gt_randomcodes_collectcodes(&fci,
                               usefirstcodes,
                               sampling_factor,
                               encseq,
                               readmode,
                               bucketkeysize,
                               (size_t)maximumspace,
                               logger,
                               timer,
                               err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (fci.differentcodes > 0 && onlyallrandomcodes)
    {
      run_allcodes_distribution(fci.allrandomcodes, fci.differentcodes);
      gt_free(fci.allrandomcodes);
      fci.allrandomcodes = NULL;
      gt_randomcodes_countocc_delete(fci.fcsl, &fci.tab);
      gt_firstcodes_spacelog_delete(fci.fcsl);
      gt_seqnumrelpos_delete(fci.buf.snrp);
      gt_Sfxmappedrangelist_delete(sfxmrlist);
      return 0;
    }
    fci.flushcount = 0;
    fci.codebuffer_total = 0;
    gt_randomcodes_fillbinsearchcache(&fci, addbscache_depth);
    if (gt_randomcodes_allocspace(&fci,
                                 numofparts,
                                 maximumspace,
                                 phase2extra,
                                 err) != 0)
    {
      haserr = true;
    }
  }
  fci.buf.nextfree = 0;
  if (!haserr)
  {
    gt_randomcodes_accumulatecounts_run(&fci,
                                       encseq,
                                       bucketkeysize,
                                       skipshorter,
                                       logger,
                                       timer);
    maxbucketsize = gt_randomcodes_partialsums(fci.fcsl, &fci.tab,
        fci.total_count);
    gt_logger_log(logger, "maximum space after computing partial sums: %.2f MB",
                  GT_MEGABYTES(gt_firstcodes_spacelog_total(fci.fcsl)));
    gt_logger_log(logger, "maxbucketsize=" GT_WU "", maxbucketsize);
    gt_randomcodes_map_sections(&fci, sfxmrlist);
    if (numofparts == 0 || maximumspace > 0)
    {
      int retval = gt_randomcodes_auto_parts(&fci,
                                            sfxmrlist,
                                            numofparts,
                                            &maximumspace,
                                            maxbucketsize,
                                            sortingdepth,
                                            totallength,
                                            maxseqlength,
                                            fci.total_count,
                                            phase2extra,
                                            err);
      if (retval < 0)
      {
        haserr = true;
      } else
      {
        numofparts = (unsigned int) retval;
      }
    }
  }
  if (!haserr)
  {
    gt_assert(numofparts > 0);
    suftabparts_rc = gt_suftabparts_rc_new(numofparts,
                                     NULL,
                                     &fci.tab,
                                     sfxmrlist,
                                     fci.total_count,
                                     0,
                                     logger);
    gt_assert(suftabparts_rc != NULL);
    gt_suftabparts_rc_showallrecords(suftabparts_rc, true);
    gt_Sfxmappedrangelist_delete(sfxmrlist);
    sfxmrlist = NULL;
    gt_randomcodes_samples_delete(fci.fcsl, &fci.tab);
    gt_assert(fci.buf.nextfree == 0);
    gt_assert(fci.mappedleftborder != NULL);
    gt_Sfxmappedrange_usetmp(fci.mappedleftborder,
        gt_randomcodes_outfilenameleftborder(&fci.tab),
        (void **)
        gt_randomcodes_leftborder_address(&fci.tab),
        gt_randomcodes_leftborder_entries(&fci.tab),
        true);
    if (gt_suftabparts_rc_numofparts(suftabparts_rc) > 1U)
    {
      gt_assert(fci.allrandomcodes != NULL);
      gt_assert(fci.mappedallrandomcodes != NULL);
      gt_Sfxmappedrange_storetmp_ulong(fci.mappedallrandomcodes,
          &fci.allrandomcodes, false);
      GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl, "allrandomcodes");
      gt_assert(fci.allrandomcodes == NULL);
    } else
    {
      gt_Sfxmappedrange_delete(fci.mappedallrandomcodes);
      fci.mappedallrandomcodes = NULL;
    }
    largest_width = gt_suftabparts_rc_largest_width(suftabparts_rc);
    fci.spmsuftab = gt_spmsuftab_new(largest_width,
                                     totallength,
                                     fci.bitsforposref,
                                     logger);
    suftab_size = gt_spmsuftab_requiredspace(largest_width, totallength,
                                             fci.bitsforposref);
    GT_FCI_ADDWORKSPACE(fci.fcsl, "suftab", suftab_size);
    fci.buf.flush_function = gt_randomcodes_insertsuffixes_flush;
    srswtab = gt_malloc(sizeof (*srswtab) * threads);
    for (threadcount = 0; threadcount < threads; threadcount++)
    {
      srswtab[threadcount]
        = gt_shortreadsort_new(maxbucketsize, maxseqlength - bucketkeysize,
                               readmode, true, false);
    }
    GT_FCI_ADDWORKSPACE(fci.fcsl, "shortreadsort",
                        threads * gt_shortreadsort_size(true, maxbucketsize,
                                                  maxseqlength-bucketkeysize));
    gt_randomcodes_allocsize_for_insertion(&fci,
                                          maximumspace,
                                          suftabparts_rc,
                                          phase2extra);
    if (!onlyaccumulation)
    {
      fci.radixsort_codepos = gt_radixsort_new_ulongpair(fci.buf.allocated);
      GT_FCI_ADDWORKSPACE(fci.fcsl, "radixsort_codepos",
                          gt_radixsort_size(fci.radixsort_codepos));
      fci.buf.spaceGtUwordPair
        = gt_radixsort_space_ulongpair(fci.radixsort_codepos);
    }
    fci.codebuffer_total = 0;
    fci.flushcount = 0;
    for (part = 0; !haserr && !onlyaccumulation &&
                   part < gt_suftabparts_rc_numofparts(suftabparts_rc); part++)
    {
      if (gt_randomcodes_process_part(&fci,
                                     encseq,
                                     readmode,
                                     bucketkeysize,
                                     sortingdepth,
                                     skipshorter,
                                     suftabparts_rc,
                                     part,
                                     maximumspace,
                                     threads,
                                     fci.total_count,
                                     withsuftabcheck,
                                     srswtab,
                                     itvprocess,
                                     itvprocess_end,
                                     itvprocessdatatab,
                                     logger,
                                     timer,
                                     err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_log_log("suftabentries=" GT_WU "", suftabentries);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "cleaning up", stdout);
  }
  if (!haserr)
  {
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl, "shortreadsort");
  }
  if (!haserr && !onlyaccumulation && srswtab != NULL)
  {
    GtUword sumofstoredvalues = 0;

    for (threadcount=0; threadcount<threads; threadcount++)
    {
      sumofstoredvalues +=
        gt_shortreadsort_sumofstoredvalues(srswtab[threadcount]);
    }
    gt_logger_log(logger, "average short read depth is %.2f",
                      (double) sumofstoredvalues/fci.total_inserted);
  }
  if (!haserr)
  {
    if (!onlyaccumulation)
    {
      gt_assert(fci.flushcount > 0);
      gt_logger_log(logger,
                    "total inserted=" GT_WU " (%.3f%% of all suffixes), "
                    "%u rounds (avg length " GT_WU ")", fci.total_inserted,
                    100.0 * (double) fci.total_inserted/totallength,
                    fci.flushcount, fci.codebuffer_total/fci.flushcount);
      gt_assert(fci.total_inserted == fci.total_count);
    } else
    {
      gt_randomcodes_delete_before_end(&fci);
    }
  }
  if (srswtab != NULL)
  {
    for (threadcount=0; threadcount<threads; threadcount++)
    {
      gt_shortreadsort_delete(srswtab[threadcount]);
    }
    gt_free(srswtab);
  }
  if (haserr)
  {
    gt_randomcodes_delete_before_end(&fci);
    gt_free(fci.allrandomcodes);
    fci.allrandomcodes = NULL;
    gt_Sfxmappedrangelist_delete(sfxmrlist);
    fci.buf.spaceGtUword = NULL;
    gt_radixsort_delete(fci.radixsort_code);
  }
  gt_suftabparts_rc_delete(suftabparts_rc);
  gt_randomcodes_countocc_delete(fci.fcsl, &fci.tab);
  gt_randomcodes_tab_delete(fci.fcsl, &fci.tab);
  if (fci.spmsuftab != NULL)
  {
    gt_spmsuftab_delete(fci.spmsuftab);
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl, "suftab");
    fci.spmsuftab = NULL;
  }
  if (fci.mappedallrandomcodes == NULL && fci.allrandomcodes != NULL)
  {
    gt_free(fci.allrandomcodes);
    fci.allrandomcodes = NULL;
    GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl, "allrandomcodes");
  }
  gt_Sfxmappedrange_delete(fci.mappedleftborder);
  gt_Sfxmappedrange_delete(fci.mappedallrandomcodes);
  gt_seqnumrelpos_delete(fci.buf.snrp);
  gt_firstcodes_spacelog_stop_diff(fci.fcsl);
  GT_FCI_SUBTRACTWORKSPACE(fci.fcsl, "encseq");
  gt_firstcodes_spacelog_delete(fci.fcsl);
  return haserr ? -1 : 0;
}
