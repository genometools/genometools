/*
  Copyright (C) 2015 Annika Seidel, <annika.seidel@studium.uni-hamburg.de>
  Copyright (C) 2015 Stefan Kurtz, <kurtz@zbh.uni-hamburg.de>
  Copyright (C) 2015 Joerg Winkler, <j.winkler@posteo.de>
  Copyright (C) 2014 Dirk Willrodt, <willrodt@zbh.uni-hamburg.de>
  Copyright (C) 2010 Sascha Steinbiss, <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2015 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>
#include "core/ma.h"
#include "core/minmax.h"
#include "core/array2dim_api.h"
#include "core/assert_api.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "match/squarededist.h"
#include "extended/alignment.h"
#include "extended/maxcoordvalue.h"
#include "extended/reconstructalignment.h"
#include "extended/squarealign.h"

#include "extended/linearalign.h"
#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

/*------------------------------global linear--------------------------------*/
static void firstEDtabRtabcolumn(GtUword *EDtabcolumn,
                                 GtUword *Rtabcolumn,
                                 GtUword ulen,
                                 GtUword gapcost)
{
  GtUword rowindex;
  EDtabcolumn[0] = 0;
  Rtabcolumn[0]  = 0;

  for (rowindex=1; rowindex <= ulen; rowindex++)
  {
    EDtabcolumn[rowindex] = EDtabcolumn[rowindex-1]+ gapcost;
    Rtabcolumn[rowindex]  = rowindex;
  }
}

static void nextEDtabRtabcolumn(GtUword *EDtabcolumn,
                                GtUword *Rtabcolumn,
                                GtUword colindex,
                                GtUword midcolumn,
                                GtUchar b,
                                const GtUchar *useq,
                                GtUword ustart,
                                GtUword ulen,
                                const GtScoreHandler *scorehandler)
{
  GtUword rowindex, val, gapcost,
          northwestEDtabentry,
          westEDtabentry,
          northwestRtabentry,
          westRtabentry = 0;

  gt_assert(scorehandler);
  gapcost = gt_scorehandler_get_gapscore(scorehandler);

  westEDtabentry = EDtabcolumn[0];
  EDtabcolumn[0] += gapcost;

  if (colindex > midcolumn)
  {
    Rtabcolumn[0] = 0;
  }

  for (rowindex = 1UL; rowindex <= ulen; rowindex++) {
    northwestEDtabentry = westEDtabentry;
    northwestRtabentry  = westRtabentry;
    westEDtabentry = EDtabcolumn[rowindex];
    westRtabentry  = Rtabcolumn[rowindex];
    EDtabcolumn[rowindex] += gapcost; /* 1. recurrence */

    /* 2. recurrence: */
    val = northwestEDtabentry + gt_scorehandler_get_replacement(scorehandler,
                                                    useq[ustart+rowindex-1], b);
    if (val <= EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (colindex > midcolumn)
      {
        Rtabcolumn[rowindex] = northwestRtabentry;
      }
    }
    /* 3. recurrence: */
    if ((val = EDtabcolumn[rowindex-1]+gapcost) < EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (colindex > midcolumn)
      {
        Rtabcolumn[rowindex] = Rtabcolumn[rowindex-1];
      }
    }
  }
}

static GtUword evaluateallEDtabRtabcolumns(GtUword *EDtabcolumn,
                                           GtUword *Rtabcolumn,
                                           const GtScoreHandler *scorehandler,
                                           GtUword midcol,
                                           const GtUchar *useq,
                                           GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vstart,
                                           GtUword vlen)
{
  GtUword gapcost, colindex;
  gt_assert(scorehandler && EDtabcolumn && Rtabcolumn);

  gapcost = gt_scorehandler_get_gapscore(scorehandler);
  firstEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, ulen, gapcost);

  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, colindex, midcol,
                        vseq[vstart+colindex-1], useq, ustart,
                        ulen, scorehandler);
  }
  return EDtabcolumn[ulen];
}

static void determineCtab0(GtUword *Ctab, const GtScoreHandler *scorehandler,
                           GtUchar vseq0, const GtUchar *useq, GtUword ustart)
{
  GtUword rowindex, repl, mincost = GT_UWORD_MAX;

  if (Ctab[1] == 0)
  {
    Ctab[0] = 0;
    return;
  }
  for (rowindex = 0; rowindex < Ctab[1]; rowindex++)
  {
    repl = gt_scorehandler_get_replacement(scorehandler,
                           vseq0, useq[ustart+rowindex]);

    if (repl == 0)
    {
      Ctab[0] = rowindex;
      return;
    }
    if (repl <= mincost)
    {
      mincost = repl;
      Ctab[0] = rowindex;
    }
  }

  if (mincost > 2 * gt_scorehandler_get_gapscore(scorehandler))
  {
    Ctab[0] = (Ctab[1] > 0) ?  Ctab[1]-1 : 0;
  }
}

#ifdef GT_THREADS_ENABLED
typedef struct{
  GtLinspaceManagement *spacemanager;
  const GtScoreHandler *scorehandler;
  const GtUchar        *useq, *vseq;
  GtUword              ustart, ulen, vstart, vlen,
                       *Ctab, rowoffset,
                       threadidx, /* ensures threads do not overlap */
                       *threadcount;
}GtLinearCrosspointthreadinfo;

static GtLinearCrosspointthreadinfo
            set_LinearCrosspointthreadinfo(GtLinspaceManagement *spacemanager,
                                           const GtScoreHandler *scorehandler,
                                           const GtUchar *useq,
                                           GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vstart,
                                           GtUword vlen,
                                           GtUword *Ctab,
                                           GtUword rowoffset,
                                           GtUword threadidx,
                                           GtUword *threadcount)
{
  GtLinearCrosspointthreadinfo threadinfo;
  threadinfo.spacemanager = spacemanager;
  threadinfo.scorehandler = scorehandler;
  threadinfo.useq = useq;
  threadinfo.ustart = ustart;
  threadinfo.ulen = ulen;
  threadinfo.vseq = vseq;
  threadinfo.vstart = vstart;
  threadinfo.vlen = vlen;
  threadinfo.Ctab = Ctab;
  threadinfo.rowoffset = rowoffset;
  threadinfo.threadidx = threadidx;
  threadinfo.threadcount = threadcount;

  return threadinfo;
}
static GtUword evaluatelinearcrosspoints(GtLinspaceManagement *spacemanager,
                                         const GtScoreHandler *scorehandler,
                                         const GtUchar *useq,
                                         GtUword ustart,
                                         GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart,
                                         GtUword vlen,
                                         GtUword *Ctab,
                                         GtUword rowoffset,
                                         GT_UNUSED GtUword threadidx,
                                         GT_UNUSED GtUword *threadcount);

static void *evaluatelinearcrosspoints_thread_caller(void *data)
{
  GtLinearCrosspointthreadinfo *threadinfo =
                                         (GtLinearCrosspointthreadinfo *) data;
  (void) evaluatelinearcrosspoints(threadinfo->spacemanager,
                                   threadinfo->scorehandler,
                                   threadinfo->useq,
                                   threadinfo->ustart,
                                   threadinfo-> ulen,
                                   threadinfo->vseq,
                                   threadinfo->vstart,
                                   threadinfo->vlen,
                                   threadinfo->Ctab,
                                   threadinfo->rowoffset,
                                   threadinfo->threadidx,
                                   threadinfo->threadcount);
  return NULL;
}
#endif

/* evaluate crosspoints in recursive way */
static GtUword evaluatelinearcrosspoints(GtLinspaceManagement *spacemanager,
                                         const GtScoreHandler *scorehandler,
                                         const GtUchar *useq,
                                         GtUword ustart, GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart, GtUword vlen,
                                         GtUword *Ctab,
                                         GtUword rowoffset,
                                         GT_UNUSED GtUword threadidx,
                                         GT_UNUSED GtUword *threadcount)
{
  GtUword midrow, midcol, distance, *EDtabcolumn = NULL, *Rtabcolumn = NULL;
#ifdef GT_THREADS_ENABLED
  GtThread *t1 = NULL, *t2 = NULL;
  GtLinearCrosspointthreadinfo threadinfo1, threadinfo2;
#endif

  if (vlen >= 2UL)
  {
    if (ulen == 0)
    {
      GtUword i;
      for (i = 0; i <= vlen; i++)
        Ctab[i] = rowoffset;
      return rowoffset;
    }

#ifdef GT_THREADS_ENABLED
    if (gt_jobs == 1)
    {
#endif
      if (gt_linspace_management_checksquare(spacemanager, ulen,vlen,
                                             sizeof (GtUword),
                                             sizeof (Rtabcolumn)))
      { /* product of subsquences is lower than space allocated already or
         * lower than timesquarfactor * ulen*/
        return gt_squarealign_ctab(spacemanager, scorehandler, Ctab, useq,
                                   ustart, ulen, vseq, vstart, vlen, rowoffset);
      }
#ifdef GT_THREADS_ENABLED
    }
#endif

    midcol = GT_DIV2(vlen);
    Rtabcolumn = gt_linspace_management_get_rTabspace(spacemanager);
    EDtabcolumn = gt_linspace_management_get_valueTabspace(spacemanager);
    Rtabcolumn = Rtabcolumn + rowoffset + threadidx;
    EDtabcolumn = EDtabcolumn + rowoffset + threadidx;

    distance = evaluateallEDtabRtabcolumns(EDtabcolumn, Rtabcolumn,
                                           scorehandler, midcol,
                                           useq, ustart, ulen,
                                           vseq, vstart, vlen);
    midrow = Rtabcolumn[ulen];
    Ctab[midcol] = rowoffset + midrow;

#ifdef GT_THREADS_ENABLED
    if (*threadcount + 2 > gt_jobs)
    {
#endif
      /* upper left corner */
      (void) evaluatelinearcrosspoints(spacemanager, scorehandler,
                                       useq, ustart, midrow,
                                       vseq, vstart, midcol,
                                       Ctab, rowoffset,
                                       threadidx, threadcount);

      /* bottom right corner */
      (void) evaluatelinearcrosspoints(spacemanager, scorehandler,
                                       useq, ustart + midrow,
                                       ulen - midrow,
                                       vseq, vstart + midcol,
                                       vlen - midcol,
                                       Ctab + midcol,
                                       rowoffset + midrow,
                                       threadidx, threadcount);
#ifdef GT_THREADS_ENABLED
    }
    else
    {
      threadinfo1 = set_LinearCrosspointthreadinfo(spacemanager, scorehandler,
                                                   useq, ustart, midrow,
                                                   vseq, vstart, midcol,
                                                   Ctab, rowoffset,
                                                   threadidx, threadcount);
      (*threadcount)++;
      t1 = gt_thread_new(evaluatelinearcrosspoints_thread_caller,
                         &threadinfo1, NULL);

      threadinfo2 = set_LinearCrosspointthreadinfo(spacemanager, scorehandler,
                                                   useq, ustart + midrow,
                                                   ulen - midrow,
                                                   vseq, vstart + midcol,
                                                   vlen - midcol,
                                                   Ctab + midcol,
                                                   rowoffset + midrow,
                                                   threadidx + GT_DIV2(midcol),
                                                   threadcount);
      (*threadcount)++;
      t2 = gt_thread_new(evaluatelinearcrosspoints_thread_caller,
                         &threadinfo2, NULL);

      gt_thread_join(t1);
      (*threadcount)--;
      gt_thread_join(t2);
      (*threadcount)--;
      gt_thread_delete(t1);
      gt_thread_delete(t2);
    }
#endif
    return distance;
  }
  return 0;
}

/* calculating alignment in linear space */
GtUword gt_calc_linearalign(GtLinspaceManagement *spacemanager,
                            const GtScoreHandler *scorehandler,
                            GtAlignment *align,
                            const GtUchar *useq,
                            GtUword ustart,
                            GtUword ulen,
                            const GtUchar *vseq,
                            GtUword vstart,
                            GtUword vlen)
{
  GtUword distance, gapcost, *Ctab, *EDtabcolumn, *Rtabcolumn, threadcount = 1;

  gt_assert(scorehandler);
  gt_linspace_management_set_ulen(spacemanager,ulen);
  gapcost = gt_scorehandler_get_gapscore(scorehandler);

  if (ulen == 0UL)
  {
    return gt_reconstructalignment_trivial_insertion(align, vlen, gapcost);
  }
  else if (vlen == 0UL)
  {
    return gt_reconstructalignment_trivial_deletion(align, ulen, gapcost);
  }
  else if (vlen == 1UL)
  {
    gt_linspace_management_check(spacemanager, (ulen+1)*(vlen+1)-1, ulen,
                                sizeof (*EDtabcolumn), sizeof (EDtabcolumn), 0);
    return gt_squarealign_calculate_generic(spacemanager, align,
                                            useq, ustart, ulen,
                                            vseq, vstart, vlen, scorehandler);
  }
  else if (gt_linspace_management_checksquare(spacemanager, ulen, vlen,
                                              sizeof (*EDtabcolumn),
                                              sizeof (*Rtabcolumn)))
  { /* call 2dim */
    return gt_squarealign_calculate_generic(spacemanager, align,
                                            useq, ustart, ulen,
                                            vseq, vstart, vlen, scorehandler);
  }

#ifdef GT_THREADS_ENABLED
  gt_linspace_management_check(spacemanager, ulen + GT_DIV2(vlen), vlen,
                               sizeof (*EDtabcolumn), sizeof (*Rtabcolumn),
                               sizeof (*Ctab));
#else
  gt_linspace_management_check(spacemanager, ulen, vlen, sizeof (*EDtabcolumn),
                               sizeof (*Rtabcolumn), sizeof (*Ctab));
#endif
  Ctab = gt_linspace_management_get_crosspointTabspace(spacemanager);

  Ctab[vlen] = ulen;
  distance = evaluatelinearcrosspoints(spacemanager, scorehandler,
                                       useq, ustart, ulen,
                                       vseq, vstart, vlen,
                                       Ctab, 0, 0, &threadcount);

  determineCtab0(Ctab, scorehandler, vseq[vstart], useq, ustart);
  gt_reconstructalignment_from_Ctab(align, Ctab, useq, ustart, vseq, vstart,
                                    vlen, scorehandler);

  return distance;
}

/* global alignment with linear gapcosts in linear space */
GtUword gt_linearalign_compute_generic(GtLinspaceManagement *spacemanager,
                                       const GtScoreHandler *scorehandler,
                                       GtAlignment *align,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen)
{
  GtUword distance;
  gt_assert(spacemanager && scorehandler && align);

  gt_alignment_set_seqs(align, useq+ustart, ulen, vseq+vstart, vlen);
  distance = gt_calc_linearalign(spacemanager, scorehandler, align,
                                 useq, ustart, ulen,
                                 vseq, vstart, vlen);

  return distance;
}

/* global alignment with linear gapcosts in linear space
   with constant cost values */
GtUword gt_linearalign_compute(GtLinspaceManagement *spacemanager,
                               GtAlignment *align,
                               const GtUchar *useq,
                               GtUword ustart,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vstart,
                               GtUword vlen,
                               GtUword matchcost,
                               GtUword mismatchcost,
                               GtUword gapcost)
{
  GtUword distance;
  GtScoreHandler *scorehandler;
  gt_assert(spacemanager && align);

  scorehandler = gt_scorehandler_new(matchcost, mismatchcost, 0, gapcost);

  distance =  gt_linearalign_compute_generic(spacemanager, scorehandler, align,
                                             useq, ustart, ulen,
                                             vseq, vstart, vlen);
  gt_scorehandler_delete(scorehandler);
  return distance;
}

/* just calculate distance, no alignment */
static void fillDPtable_unitcost(bool downcase,GtUword *dpcolumn,
                                 const GtUchar *u, GtUword ulen,
                                 const GtUchar *v, GtUword vlen)
{
  GtUword i, j , nw, we;

  for (i = 0; i <= ulen; i++)
  {
    dpcolumn[i] = i;
  }
  for (j = 1UL; j <= vlen; j++) {
    GtUchar b = downcase ? tolower((int) v[j-1]) : v[j-1];
    nw = dpcolumn[0];
    dpcolumn[0] = j;
    for (i = 1UL; i <= ulen; i++) {
      GtUchar a = downcase ? tolower((int) u[i-1]) : u[i-1];
      we = dpcolumn[i];
       /* replacement */
      dpcolumn[i] = (a == b) ? nw : (nw + 1);
      if (dpcolumn[i-1] + 1 < dpcolumn[i]) /* deletion */
      {
        dpcolumn[i] = dpcolumn[i-1] + 1;
      }
      if (we + 1 < dpcolumn[i]) /* insertion */
      {
        dpcolumn[i] = we + 1;
      }
      nw = we;
    }
  }
}

static GtUword gt_calc_linearedist(bool downcase,
                                   const GtUchar *useq, GtUword ulen,
                                   const GtUchar *vseq, GtUword vlen)
{
  GtUword *dpcolumn, edist;

  printf("%s(%*.*s,%*.*s)\n",__func__,(int) ulen,(int) ulen,(char *) useq,
                             (int) vlen,(int) vlen,(char *) vseq);
  dpcolumn = gt_malloc(sizeof *dpcolumn * (MIN(ulen,vlen) + 1));
  fillDPtable_unitcost(downcase,dpcolumn,
                       ulen <= vlen ? useq : vseq, MIN(ulen,vlen),
                       ulen <= vlen ? vseq : useq, MAX(ulen,vlen));
  edist = dpcolumn[MIN(ulen,vlen)];
  gt_free(dpcolumn);
  return edist;
}

/*-------------------------------local linear---------------------------------*/

static void firstLStabcolumn(GtWord *Ltabcolumn,
                             GtUwordPair *Starttabcolumn,
                             GtUword ulen)
{
  GtUword rowindex;

  gt_assert(Ltabcolumn != NULL && Starttabcolumn != NULL);
  for (rowindex = 0; rowindex <= ulen; rowindex++)
  {
    Ltabcolumn[rowindex] = 0;
    Starttabcolumn[rowindex].a = rowindex;
    Starttabcolumn[rowindex].b = 0;
  }
}

static void nextLStabcolumn(GtWord *Ltabcolumn,
                            GtUwordPair *Starttabcolumn,
                            const GtScoreHandler *scorehandler,
                            const GtUchar *useq,
                            GtUword ustart, GtUword ulen,
                            const GtUchar b, GtUword colindex,
                            GtMaxcoordvalue *max)
{
  GtUword rowindex;
  GtUwordPair northwestStarttabentry, westStarttabentry;
  GtWord northwestLtabentry, westLtabentry, gapscore;

  gt_assert(max != NULL);
  gapscore = gt_scorehandler_get_gapscore(scorehandler);

  westLtabentry = Ltabcolumn[0];
  westStarttabentry = Starttabcolumn[0];

  Ltabcolumn[0] = 0;
  Starttabcolumn[0].a = 0;
  Starttabcolumn[0].b = colindex;
  for (rowindex = 1UL; rowindex <= ulen; rowindex++)
  {
    GtWord val;

    northwestLtabentry = westLtabentry;
    northwestStarttabentry = westStarttabentry;
    westLtabentry = Ltabcolumn[rowindex];
    westStarttabentry = Starttabcolumn[rowindex];
    Ltabcolumn[rowindex] += gapscore;

    val = northwestLtabentry + gt_scorehandler_get_replacement(scorehandler,
                                                  useq[ustart + rowindex-1], b);

    if (val >= Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = val;
      Starttabcolumn[rowindex] = northwestStarttabentry;
    }
    if ((val = Ltabcolumn[rowindex-1]+gapscore) > Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = val;
      Starttabcolumn[rowindex]=Starttabcolumn[rowindex-1];
    }
    if (0 > Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = 0;
      Starttabcolumn[rowindex].a = rowindex;
      Starttabcolumn[rowindex].b = colindex;
    }
    if (Ltabcolumn[rowindex] > gt_maxcoordvalue_get_value(max))
    {
      gt_maxcoordvalue_coord_update(max, Ltabcolumn[rowindex],
                                    Starttabcolumn[rowindex],
                                    rowindex, colindex);
    }
  }
}

static GtMaxcoordvalue *evaluateallLScolumns(GtLinspaceManagement *spacemanager,
                                             const GtScoreHandler *scorehandler,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen)
{
  GtUword colindex;
  GtWord *Ltabcolumn;
  GtUwordPair *Starttabcolumn;
  GtMaxcoordvalue *max;

  Ltabcolumn = gt_linspace_management_get_valueTabspace(spacemanager);
  Starttabcolumn = gt_linspace_management_get_rTabspace(spacemanager);

  firstLStabcolumn(Ltabcolumn, Starttabcolumn, ulen);

  max = gt_linspace_management_get_maxspace(spacemanager);
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextLStabcolumn(Ltabcolumn, Starttabcolumn, scorehandler,
                    useq, ustart, ulen, vseq[vstart+colindex-1], colindex, max);
  }
  return max;
}

/* determining start and end of local alignment and call global function */
GtWord gt_linearalign_compute_local_generic(GtLinspaceManagement *spacemanager,
                                            const GtScoreHandler *scorehandler,
                                            GtAlignment *align,
                                            const GtUchar *useq,
                                            GtUword ustart,
                                            GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart,
                                            GtUword vlen)
{
  GtWord *Ltabcolumn,
         score = GT_WORD_MAX;
  GtUwordPair *Starttabcolumn;
  GtUword ulen_part, ustart_part, vlen_part, vstart_part;
  GtMaxcoordvalue *max;

  gt_assert(spacemanager && scorehandler && align);
  gt_linspace_management_set_ulen(spacemanager,ulen);

  if (ulen == 0UL || vlen == 0UL)
  {
     /* empty alignment */
    return 0;
  }
  else if (vlen == 1UL)
  {
    gt_linspace_management_check_local(spacemanager,
                                       (ulen+1)*(vlen+1)-1, ulen,
                                       sizeof (GtWord),
                                       sizeof (GtWord *));
    return gt_squarealign_calculate_local_generic(spacemanager, align,
                                                  useq, ustart, ulen,
                                                  vseq, vstart, vlen,
                                                  scorehandler);
  }
  else if (gt_linspace_management_checksquare_local(spacemanager, ulen, vlen,
                                                    sizeof (*Ltabcolumn),
                                                    sizeof (*Starttabcolumn)))
  {
    /* call 2dim */
    return gt_squarealign_calculate_local_generic(spacemanager, align,
                                                  useq, ustart, ulen,
                                                  vseq, vstart, vlen,
                                                  scorehandler);
  }

  gt_linspace_management_check_local(spacemanager,
                                     ulen, vlen,
                                     sizeof (*Ltabcolumn),
                                     sizeof (*Starttabcolumn));

  max = evaluateallLScolumns(spacemanager, scorehandler,
                             useq, ustart, ulen,
                             vseq, vstart, vlen);

  if (gt_maxcoordvalue_get_length_safe(max))
  {
    GtScoreHandler *costhandler;

    ustart_part = ustart+(gt_maxcoordvalue_get_start(max)).a;
    vstart_part = vstart+(gt_maxcoordvalue_get_start(max)).b;
    ulen_part = gt_maxcoordvalue_get_row_length(max);
    vlen_part = gt_maxcoordvalue_get_col_length(max);
    score = gt_maxcoordvalue_get_value(max);

    gt_alignment_set_seqs(align, useq + ustart_part, ulen_part,
                                 vseq + vstart_part, vlen_part);
    costhandler = gt_scorehandler2costhandler(scorehandler);
    /* call global function */
    gt_calc_linearalign(spacemanager,
                        costhandler, align,
                        useq, ustart_part, ulen_part,
                        vseq, vstart_part, vlen_part);
    gt_scorehandler_delete(costhandler);
  } else
  {
    /*empty alignment */
    return 0;
  }
  return score;
}

GtWord gt_linearalign_compute_local(GtLinspaceManagement *spacemanager,
                                    GtAlignment *align,
                                    const GtUchar *useq,
                                    GtUword ustart,
                                    GtUword ulen,
                                    const GtUchar *vseq,
                                    GtUword vstart,
                                    GtUword vlen,
                                    GtWord matchscore,
                                    GtWord mismatchscore,
                                    GtWord gapscore)
{
  GtWord score;
  gt_assert(align && spacemanager);
  GtScoreHandler *scorehandler = gt_scorehandler_new(matchscore,
                                                     mismatchscore,
                                                     0, gapscore);

  score = gt_linearalign_compute_local_generic(spacemanager, scorehandler,
                                               align, useq, ustart, ulen,
                                               vseq, vstart, vlen);
  gt_scorehandler_delete(scorehandler);
  return score;
}

/*-----------------------------checkfunctions--------------------------------*/

void gt_linearalign_check(GT_UNUSED bool forward,
                          const GtUchar *useq,
                          GtUword ulen,
                          const GtUchar *vseq,
                          GtUword vlen)
{
  GtAlignment *align;
  GtUword edist1, edist2, edist3, edist4,
          matchcost = 0, mismatchcost = 1, gapcost = 1;
  GtLinspaceManagement *spacemanager;
  GtScoreHandler *scorehandler;
  const bool downcase = true;

  if (memchr(useq, LINEAR_EDIST_GAP,ulen) != NULL)
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (memchr(vseq, LINEAR_EDIST_GAP,vlen) != NULL)
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  scorehandler = gt_scorehandler_new(matchcost,  mismatchcost, 0, gapcost);
  gt_scorehandler_plain(scorehandler);
  gt_scorehandler_downcase(scorehandler);
  spacemanager = gt_linspace_management_new();
  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);
  edist1 = gt_calc_linearalign(spacemanager, scorehandler, align,
                               useq, 0, ulen,
                               vseq, 0, vlen);
  edist2 = gt_squarealign_global_distance_only(useq, 0, ulen, vseq, 0, vlen,
                                               scorehandler);

  if (edist1 != edist2)
  {
    fprintf(stderr,"gt_calc_linearalign = "GT_WU" != "GT_WU
            " = gt_squarealign_global_distance_only\n", edist1,edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  edist3 = gt_alignment_eval_with_score(align, true, matchcost,
                                        mismatchcost, gapcost);

  if (edist2 != edist3)
  {
    fprintf(stderr,"gt_squarealign_global_distance_only = "GT_WU" != "GT_WU
            " = gt_alignment_eval_with_score\n", edist2,edist3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  edist4 = gt_calc_linearedist(downcase,useq, ulen, vseq, vlen);
  if (edist3 != edist4)
  {
    fprintf(stderr,"gt_alignment_eval_with_score = "GT_WU" != "GT_WU
            " = gt_calc_linearedist\n", edist3, edist4);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_linspace_management_delete(spacemanager);
  gt_scorehandler_delete(scorehandler);
  gt_alignment_delete(align);
}

void gt_linearalign_check_local(GT_UNUSED bool forward,
                                const GtUchar *useq, GtUword ulen,
                                const GtUchar *vseq, GtUword vlen)
{
  GtAlignment *align;
  GtWord score1, score2, score3, score4,
         matchscore = 2, mismatchscore = -2, gapscore = -1;
  GtLinspaceManagement *spacemanager;
  GtScoreHandler *scorehandler;

  if (memchr(useq, LINEAR_EDIST_GAP,ulen) != NULL)
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (memchr(vseq, LINEAR_EDIST_GAP,vlen) != NULL)
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  scorehandler = gt_scorehandler_new(matchscore, mismatchscore, 0, gapscore);
  gt_scorehandler_plain(scorehandler);
  spacemanager = gt_linspace_management_new();
  align = gt_alignment_new();
  score1 = gt_linearalign_compute_local_generic(spacemanager, scorehandler,
                                                align, useq, 0, ulen,
                                                vseq, 0, vlen);

  score2 = gt_alignment_eval_with_score(align, true, matchscore,
                                        mismatchscore, gapscore);
  gt_linspace_management_delete(spacemanager);
  gt_scorehandler_delete(scorehandler);
  if (score1 != score2)
  {
    fprintf(stderr,"gt_linearalign_compute_local_generic = "GT_WD" != "GT_WD
            " = gt_alignment_eval_generic_with_score\n", score1, score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_alignment_reset(align);
  score3 = gt_squarealign_calculate_local(NULL, align, useq, 0, ulen,
                                          vseq, 0, vlen, matchscore,
                                          mismatchscore, gapscore);

  if (score1 != score3)
  {
    fprintf(stderr,"gt_linearalign_compute_local_generic = "GT_WD" != "GT_WD
            " = gt_squarealign_calculate_local\n", score1, score3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  score4 = gt_alignment_eval_with_score(align, true, matchscore,
                                        mismatchscore, gapscore);
  if (score3 != score4)
  {
    fprintf(stderr,"gt_squarealign_calculate_local = "GT_WD" != "GT_WD
            " = gt_alignment_eval_generic_with_score\n", score3, score4);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_alignment_delete(align);
}
