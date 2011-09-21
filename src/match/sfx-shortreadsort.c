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

#include <stdint.h>
#include "core/intbits.h"
#include "core/stack-inlined.h"
#include "core/minmax.h"
#include "core/assert_api.h"
#include "sfx-lcpvalues.h"
#include "sfx-shortreadsort.h"
#include "spmsuftab.h"

#ifdef _LP64
#define GT_NUMOFTBEVALUEFOR100 3
#else
#define GT_NUMOFTBEVALUEFOR100 6
#endif

typedef struct
{
  unsigned long suffix;
  unsigned long seqnum_relpos;
  GtTwobitencoding tbe[GT_NUMOFTBEVALUEFOR100];
  unsigned int unitsnotspecial;
} GtShortreadsort;

typedef uint32_t GtShortreadsortreftype;

struct GtShortreadsortworkinfo
{
  GtShortreadsortreftype *shortreadsortrefs;
  GtShortreadsort *shortreadsortinfo;
  GtLcpvalues *tableoflcpvalues; /* always NULL in the context of
                                    firstcodes; otherwise: NULL iff
                                    lcpvalues are not required. */
  uint16_t *firstcodeslcpvalues;
  unsigned long numofentries,
                tmplcplen;
  unsigned int bitsforrelpos;
  bool fwd, complement;
};

size_t gt_shortreadsort_size(unsigned long maxvalue)
{
  return (sizeof (GtShortreadsort) +
         sizeof (GtShortreadsortreftype)) * (maxvalue+1);
}

const uint16_t *gt_shortreadsort_lcpvalues(const GtShortreadsortworkinfo *srsw)
{
  gt_assert(srsw->firstcodeslcpvalues != NULL);

  return srsw->firstcodeslcpvalues;
}

GtShortreadsortworkinfo *gt_shortreadsort_new(unsigned long maxwidth,
                                              GtReadmode readmode,
                                              unsigned int bitsforrelpos,
                                              bool firstcodes)
{
  unsigned long idx;
  GtShortreadsortworkinfo *srsw;
#ifndef NDEBUG
  const unsigned long maxshortreadsort
   = sizeof (GtShortreadsortreftype) == sizeof (uint16_t) ? UINT16_MAX
                                                          : UINT32_MAX;
#endif

  srsw = gt_malloc(sizeof(*srsw));
  gt_assert(maxwidth <= maxshortreadsort);
  srsw->numofentries = maxwidth + 1;
  srsw->bitsforrelpos = bitsforrelpos;
  srsw->shortreadsortinfo
    = gt_malloc(sizeof (*srsw->shortreadsortinfo) * srsw->numofentries);
  srsw->shortreadsortrefs
    = gt_malloc(sizeof (*srsw->shortreadsortrefs) * srsw->numofentries);
  if (firstcodes)
  {
    srsw->firstcodeslcpvalues
      = gt_malloc(sizeof (*srsw->firstcodeslcpvalues) * srsw->numofentries);
    srsw->firstcodeslcpvalues[0] = 0; /* since it is not set otherwise */
  } else
  {
    srsw->firstcodeslcpvalues = NULL;
  }
  srsw->fwd = GT_ISDIRREVERSE(readmode) ? false : true;
  srsw->complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;
  srsw->tableoflcpvalues = NULL;
  for (idx = 0; idx < srsw->numofentries; idx++)
  {
    srsw->shortreadsortrefs[idx] = (GtShortreadsortreftype) idx;
  }
  return srsw;
}

void gt_shortreadsort_delete(GtShortreadsortworkinfo *srsw)
{
  if (srsw != NULL)
  {
    gt_free(srsw->shortreadsortinfo);
    srsw->shortreadsortinfo = NULL;
    gt_free(srsw->shortreadsortrefs);
    srsw->shortreadsortrefs = NULL;
    gt_free(srsw->firstcodeslcpvalues);
    srsw->firstcodeslcpvalues = NULL;
    gt_free(srsw);
  }
}

void gt_shortreadsort_assigntableoflcpvalues(
          GtShortreadsortworkinfo *srsw,GtLcpvalues *tableoflcpvalues)
{
  srsw->tableoflcpvalues = tableoflcpvalues;
}

static int gt_shortreadsort_compare(const GtShortreadsort *aq,
                                    const GtShortreadsort *bq,
                                    GtShortreadsortworkinfo *srsw)
{
  int idx, retval;
  unsigned int maxprefix;
  GtCommonunits commonunits;

  for (idx=0, maxprefix = (unsigned int) GT_UNITSIN2BITENC;
       idx<GT_NUMOFTBEVALUEFOR100;
       idx++, maxprefix+=(unsigned int) GT_UNITSIN2BITENC)
  {
    if (aq->unitsnotspecial >= maxprefix &&
        bq->unitsnotspecial >= maxprefix)
    {
      if (aq->tbe[idx] != bq->tbe[idx])
      {
        retval = gt_encseq_compare_pairof_different_twobitencodings(
                                                     srsw->fwd,
                                                     srsw->complement,
                                                     &commonunits,
                                                     aq->tbe[idx],
                                                     bq->tbe[idx]);
        srsw->tmplcplen
          = (unsigned long) (maxprefix - GT_UNITSIN2BITENC +
                             commonunits.common);
        return retval;
      }
    } else
    {
      GtEndofTwobitencoding tbe_a, tbe_b;

      tbe_a.tbe = aq->tbe[idx];
      tbe_b.tbe = bq->tbe[idx];
      tbe_a.position = aq->suffix;
      tbe_b.position = bq->suffix;
      tbe_a.unitsnotspecial
        = aq->unitsnotspecial >= maxprefix
           ? maxprefix
           : aq->unitsnotspecial + GT_UNITSIN2BITENC - maxprefix;
      tbe_b.unitsnotspecial
        = bq->unitsnotspecial >= maxprefix
           ? maxprefix
           : bq->unitsnotspecial + GT_UNITSIN2BITENC - maxprefix;
      retval = gt_encseq_compare_pairof_twobitencodings(srsw->fwd,
                                                        srsw->complement,
                                                        &commonunits,
                                                        &tbe_a,
                                                        &tbe_b);
      srsw->tmplcplen = (unsigned long) (maxprefix - GT_UNITSIN2BITENC +
                                         commonunits.common);
      return retval;
    }
  }
  gt_assert(false);
  return 0;
}

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME) shortread_##NAME

#define shortread_ARRAY_GET(ARR,IDX)\
        (unsigned long) data->shortreadsortrefs[IDX]

#define shortread_ARRAY_SET(ARR,IDX,VALUE)\
        data->shortreadsortrefs[IDX] = (GtShortreadsortreftype) VALUE

typedef GtShortreadsortworkinfo * QSORTNAME(Datatype);

static int QSORTNAME(qsortcmparr) (unsigned long a,
                                   unsigned long b,
                                   const QSORTNAME(Datatype) data)
{
  return gt_shortreadsort_compare(
                      data->shortreadsortinfo + QSORTNAME(ARRAY_GET)(NULL,a),
                      data->shortreadsortinfo + QSORTNAME(ARRAY_GET)(NULL,b),
                      data);
}

typedef unsigned long QSORTNAME(Sorttype);

/*
 * Qsort routine from Bentley & McIlroy's ``Engineering a Sort Function''.
 */

#ifndef GT_QSORT_ARR_SWAP
#define GT_QSORT_ARR_SWAP(ARR,A,B)\
        if ((A) != (B))\
        {\
          tmp = QSORTNAME(ARRAY_GET)(ARR,A);\
          QSORTNAME(ARRAY_SET)(ARR,A,QSORTNAME(ARRAY_GET)(ARR,B));\
          QSORTNAME(ARRAY_SET)(ARR,B,tmp);\
        }
#endif

#ifndef GT_QSORT_ARR_VECSWAP
#define GT_QSORT_ARR_VECSWAP(ARR,A,B,N)\
        aidx = A;\
        bidx = B;\
        while ((N)-- > 0)\
        {\
          tmp = QSORTNAME(ARRAY_GET)(ARR,aidx);\
          QSORTNAME(ARRAY_SET)(ARR,aidx,QSORTNAME(ARRAY_GET)(ARR,bidx));\
          QSORTNAME(ARRAY_SET)(ARR,bidx,tmp);\
          aidx++;\
          bidx++;\
        }
#endif

static inline unsigned long QSORTNAME(gt_inlined_qsort_arr_r_med3)
                     (unsigned long a, unsigned long b, unsigned long c,
                      QSORTNAME(Datatype) data)
{
  return QSORTNAME(qsortcmparr) (a, b, data) < 0
           ? (QSORTNAME(qsortcmparr) (b, c, data) < 0
                ? b
                : (QSORTNAME(qsortcmparr) (a, c, data) < 0
                     ? c : a))
           : (QSORTNAME(qsortcmparr) (b, c, data) > 0
                ? b
                : (QSORTNAME(qsortcmparr) (a, c, data) < 0
                     ? a
                     : c));
}

#ifndef GT_STACK_INTERVALARRAYTOBESORTED_DEFINED
typedef struct
{
  unsigned long startindex,
                len;
} Intervalarrtobesorted;

GT_STACK_DECLARESTRUCT(Intervalarrtobesorted,32UL);
#define GT_STACK_INTERVALARRAYTOBESORTED_DEFINED
#endif

static void QSORTNAME(gt_inlinedarr_qsort_r) (
                                   unsigned long insertionsortthreshold,
                                   bool handlenotswapped,
                                   unsigned long len,
                                   QSORTNAME(Datatype) data,
                                   unsigned long depth,
                                   unsigned long subbucketleft)
{
  unsigned long tmp, pa, pb, pc, pd, pl, pm, pn, aidx, bidx, s,
                smallermaxlcp, greatermaxlcp;
  int r;
  bool swapped;
  GtStackIntervalarrtobesorted intervalstack;
  Intervalarrtobesorted current;

  GT_STACK_INIT(&intervalstack,32UL);
  current.startindex = 0;
  current.len = len;
  GT_STACK_PUSH(&intervalstack,current);
  if (insertionsortthreshold <= 2UL)
  {
    insertionsortthreshold = 6UL;
  }
  while (!GT_STACK_ISEMPTY(&intervalstack))
  {
    swapped = false;
    current = GT_STACK_POP(&intervalstack);
    if (current.len <= insertionsortthreshold)
    {
      for (pm = current.startindex + 1;
           pm < current.startindex + current.len; pm++)
      {
        for (pl = pm; pl > current.startindex; pl--)
        {
          r = QSORTNAME(qsortcmparr) (pl - 1, pl, data);
          if (data->tableoflcpvalues != NULL)
          {
            unsigned long lcpindex = subbucketleft + pl;
            if (pl < pm && r > 0)
            {
              lcptab_update(data->tableoflcpvalues,lcpindex+1,
                            lcpsubtab_getvalue(data->tableoflcpvalues,
                                               lcpindex));
            }
            lcptab_update(data->tableoflcpvalues,lcpindex,
                          depth + data->tmplcplen);
          } else
          {
            if (data->firstcodeslcpvalues != NULL)
            {
              if (pl < pm && r > 0)
              {
                data->firstcodeslcpvalues[pl+1] = data->firstcodeslcpvalues[pl];
              }
              gt_assert(depth + data->tmplcplen <= UINT16_MAX);
              data->firstcodeslcpvalues[pl]
                = (uint16_t) (depth + data->tmplcplen);
            }
          }
          if (r <= 0)
          {
            break;
          }
          GT_QSORT_ARR_SWAP (arr, pl, pl - 1);
        }
      }
      continue;
    }
    pm = current.startindex + GT_DIV2 (current.len);
    if (current.len > 7UL)
    {
      pl = current.startindex;
      pn = current.startindex + current.len - 1;
      if (current.len > 40UL)
      {
        s = GT_DIV8 (current.len);
        pl = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pl, pl + s,
                                                     pl + GT_MULT2 (s), data);
        gt_assert(pm >= s);
        pm = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pm - s, pm,
                                                     pm + s, data);
        gt_assert(pn >= GT_MULT2(s));
        pn = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pn - GT_MULT2 (s),
                                                     pn - s, pn, data);
      }
      pm = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pl, pm, pn, data);
    }
    GT_QSORT_ARR_SWAP (arr, current.startindex, pm);
    pa = pb = current.startindex + 1;
    pc = pd = current.startindex + current.len - 1;
    smallermaxlcp = greatermaxlcp = 0;
    while (1)
    {
      while (pb <= pc)
      {
        r = QSORTNAME(qsortcmparr) (pb, current.startindex, data);
        if (r > 0)
        {
          GT_BSR_UPDATEMAXLCP(greatermaxlcp,data->tmplcplen);
          break;
        }
        if (r == 0)
        {
          swapped = true;
          GT_QSORT_ARR_SWAP (arr, pa, pb);
          pa++;
        } else
        {
          GT_BSR_UPDATEMAXLCP(smallermaxlcp,data->tmplcplen);
        }
        pb++;
      }
      while (pb <= pc)
      {
        r = QSORTNAME(qsortcmparr) (pc, current.startindex, data);
        if (r < 0)
        {
          GT_BSR_UPDATEMAXLCP(smallermaxlcp,data->tmplcplen);
          break;
        }
        if (r == 0)
        {
          swapped = true;
          GT_QSORT_ARR_SWAP (arr, pc, pd);
          gt_assert(pd > 0);
          pd--;
        } else
        {
          GT_BSR_UPDATEMAXLCP(greatermaxlcp,data->tmplcplen);
        }
        gt_assert(pc > 0);
        pc--;
      }
      if (pb > pc)
      {
        break;
      }
      GT_QSORT_ARR_SWAP (arr, pb, pc);
      swapped = true;
      pb++;
      gt_assert(pc > 0);
      pc--;
    }
    /* The following switch is not explained in the above mentioned
       paper and therefore we ignore it. */
    if (handlenotswapped && !swapped)
    {                                  /* Switch to insertion sort */
      gt_assert(current.len <= 40UL);
      for (pm = current.startindex + 1;
           pm < current.startindex + current.len; pm++)
      {
        for (pl = pm; pl > current.startindex; pl--)
        {
          r = QSORTNAME(qsortcmparr) (pl - 1, pl, data);
          if (r <= 0)
          {
            break;
          }
          GT_QSORT_ARR_SWAP (arr, pl, pl - 1);
        }
      }
      continue;
    }
    pn = current.startindex + current.len;
    gt_assert(pa >= current.startindex);
    gt_assert(pb >= pa);
    s = MIN ((unsigned long) (pa - current.startindex),
             (unsigned long) (pb - pa));
    gt_assert(pb >= s);
    GT_QSORT_ARR_VECSWAP (arr, current.startindex, pb - s, s);
    gt_assert(pd >= pc);
    gt_assert(pn > pd);
    s = MIN ((unsigned long) (pd - pc), (unsigned long) (pn - pd - 1));
    gt_assert(pn > s);
    GT_QSORT_ARR_VECSWAP (arr, pb, pn - s, s);
    gt_assert(pb >= pa);
    if ((s = (unsigned long) (pb - pa)) > 0)
    {
      if (data->tableoflcpvalues != NULL)
      {
        /*
          left part has suffix with lcp up to length smallermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the left
          which is at a minimum distance to the pivot and thus to an
          element in the final part of the left side.
        */
        lcptab_update(data->tableoflcpvalues,subbucketleft+current.startindex+s,
                      depth+smallermaxlcp);
      } else
      {
        if (data->firstcodeslcpvalues != NULL)
        {
          gt_assert(depth + smallermaxlcp <= UINT16_MAX);
          data->firstcodeslcpvalues[current.startindex+s]
            = (uint16_t) (depth + smallermaxlcp);
        }
      }
      if (s > 1UL)
      {
        current.len = s;
        GT_STACK_PUSH(&intervalstack,current);
      }
    }
    gt_assert(pd >= pc);
    if ((s = (unsigned long) (pd - pc)) > 0)
    {
      if (data->tableoflcpvalues != NULL)
      {
        /*
          right part has suffix with lcp up to length largermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the right
          which is at a minimum distance to the pivot and thus to an
          element in the first part of the right side.
        */
        gt_assert(pn >= s);
        lcptab_update(data->tableoflcpvalues, subbucketleft + pn - s,
                      depth + greatermaxlcp);
      } else
      {
        if (data->firstcodeslcpvalues != NULL)
        {
          gt_assert(depth + greatermaxlcp <= UINT16_MAX);
          data->firstcodeslcpvalues[pn - s]
            = (uint16_t) (depth + greatermaxlcp);
        }
      }
      if (s > 1UL)
      {
        gt_assert(pn >= s);
        current.startindex = pn - s;
        current.len = s;
        GT_STACK_PUSH(&intervalstack,current);
      }
    }
  }
  GT_STACK_DELETE(&intervalstack);
}

void gt_shortreadsort_sssp_sort(GtShortreadsortworkinfo *srsw,
                                const GtEncseq *encseq,
                                GtReadmode readmode,
                                GtEncseqReader *esr,
                                GtSuffixsortspace *sssp,
                                unsigned long subbucketleft,
                                unsigned long width,
                                unsigned long depth)
{
  unsigned long idx, pos;
  GtSuffixsortspace_exportptr *exportptr;

  gt_assert(srsw->bitsforrelpos == 0); /* this is only needed in array_sort */
  exportptr = gt_suffixsortspace_exportptr(subbucketleft, sssp);
  if (exportptr->ulongtabsectionptr != NULL)
  {
    for (idx = 0; idx < width; idx++)
    {
      pos = exportptr->ulongtabsectionptr[idx];
      srsw->shortreadsortinfo[idx].suffix = pos;
      srsw->shortreadsortinfo[idx].unitsnotspecial
        = gt_encseq_extract2bitencvector(srsw->shortreadsortinfo[idx].tbe,
                                         GT_NUMOFTBEVALUEFOR100,
                                         encseq,
                                         esr,
                                         readmode,
                                         pos+depth);
    }
  } else
  {
    for (idx = 0; idx < width; idx++)
    {
      pos = (unsigned long) exportptr->uinttabsectionptr[idx];
      srsw->shortreadsortinfo[idx].suffix = pos;
      srsw->shortreadsortinfo[idx].unitsnotspecial
        = gt_encseq_extract2bitencvector(srsw->shortreadsortinfo[idx].tbe,
                                         GT_NUMOFTBEVALUEFOR100,
                                         encseq,
                                         esr,
                                         readmode,
                                         pos+depth);
    }
  }
  QSORTNAME(gt_inlinedarr_qsort_r) (6UL, false, width, srsw, depth,
                                    subbucketleft);
  if (exportptr->ulongtabsectionptr != NULL)
  {
    for (idx = 0; idx < width; idx++)
    {
      exportptr->ulongtabsectionptr[idx]
        = srsw->shortreadsortinfo[srsw->shortreadsortrefs[idx]].suffix;
      srsw->shortreadsortrefs[idx] = (GtShortreadsortreftype) idx;
      if (exportptr->ulongtabsectionptr[idx] == 0)
      {
        gt_suffixsortspace_updatelongest(sssp,idx);
      }
    }
  } else
  {
    for (idx = 0; idx < width; idx++)
    {
      exportptr->uinttabsectionptr[idx]
        = (uint32_t) srsw->shortreadsortinfo[srsw->shortreadsortrefs[idx]]
                                            .suffix;
      srsw->shortreadsortrefs[idx] = (GtShortreadsortreftype) idx;
      if (exportptr->uinttabsectionptr[idx] == 0)
      {
        gt_suffixsortspace_updatelongest(sssp,idx);
      }
    }
  }
  gt_suffixsortspace_export_done(sssp);
}

#define GT_COMBINE_SEQNUM_RELPOS(SEQNUM,SHIFT,RELPOS)\
        ((SEQNUM) << (SHIFT) | (RELPOS))

void gt_shortreadsort_array_sort(unsigned long *suftab_bucket,
                                 unsigned long *seqnum_relpos_bucket,
                                 GtShortreadsortworkinfo *srsw,
                                 const GtEncseq *encseq,
                                 GtSpmsuftab *spmsuftab,
                                 unsigned long subbucketleft,
                                 unsigned long width,
                                 unsigned long depth)
{
  unsigned long idx, pos, seqnum, relpos;

  for (idx = 0; idx < width; idx++)
  {
    srsw->shortreadsortrefs[idx] = (GtShortreadsortreftype) idx;
  }
  for (idx = 0; idx < width; idx++)
  {
    pos = gt_spmsuftab_get(spmsuftab,subbucketleft + idx);
    seqnum = gt_encseq_seqnum(encseq,pos+depth);
    relpos = pos + depth - gt_encseq_seqstartpos(encseq,seqnum);
    srsw->shortreadsortinfo[idx].suffix = pos;
    srsw->shortreadsortinfo[idx].seqnum_relpos
      = GT_COMBINE_SEQNUM_RELPOS(seqnum,srsw->bitsforrelpos,relpos);
    srsw->shortreadsortinfo[idx].unitsnotspecial
      = gt_encseq_relpos_extract2bitencvector(
                         srsw->shortreadsortinfo[idx].tbe,
                         GT_NUMOFTBEVALUEFOR100,
                         encseq,
                         seqnum,
                         relpos);
  }
  QSORTNAME(gt_inlinedarr_qsort_r) (6UL, false, width, srsw, depth,
                                    subbucketleft);
  for (idx = 0; idx < width; idx++)
  {
    suftab_bucket[idx]
      = srsw->shortreadsortinfo[srsw->shortreadsortrefs[idx]].suffix;
    seqnum_relpos_bucket[idx]
      = srsw->shortreadsortinfo[srsw->shortreadsortrefs[idx]].seqnum_relpos;
  }
}
