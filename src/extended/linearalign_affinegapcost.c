/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/error.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif
#include "core/types_api.h"
#include "extended/affinealign.h"
#include "extended/maxcoordvalue.h"
#include "extended/reconstructalignment.h"

#include "extended/linearalign_affinegapcost.h"
#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

typedef struct {
  GtUwordPair Rstart, Dstart, Istart;
} Starttabentry;

/*-------------------------------global affine--------------------------------*/
GtAffineAlignEdge gt_linearalign_affinegapcost_set_edge(GtWord Rdist,
                                                        GtWord Ddist,
                                                        GtWord Idist)
{
  GtUword minvalue = MIN3(Rdist, Ddist, Idist);

  if (Rdist == minvalue)
    return Affine_R;
  else if (Ddist == minvalue)
    return Affine_D;
  else if (Idist == minvalue)
    return Affine_I;

  return Affine_X;
}

static inline GtAffineAlignRnode get_Rtabentry(const GtAffineAlignRtabentry
                                                *rtab, GtAffineAlignEdge edge)
{
  if (edge == Affine_R)
  {
    return rtab->val_R;
  }
  if (edge == Affine_D)
  {
    return rtab->val_D;
  }
  gt_assert(edge == Affine_I);
  return rtab->val_I;
}

static inline void firstAtabRtabentry(GtAffinealignDPentry *Atabcolumn,
                                      GtUword gap_opening,
                                      GtAffineAlignEdge edge)
{
  Atabcolumn[0].Redge = Affine_X;
  Atabcolumn[0].Dedge = Affine_X;
  Atabcolumn[0].Iedge = Affine_X;

  switch (edge) {
  case Affine_R:
    Atabcolumn[0].Rvalue = 0;
    Atabcolumn[0].Dvalue = GT_WORD_MAX;
    Atabcolumn[0].Ivalue = GT_WORD_MAX;
    break;
  case Affine_D:
    Atabcolumn[0].Rvalue = GT_WORD_MAX;
    Atabcolumn[0].Dvalue = 0;
    Atabcolumn[0].Ivalue = GT_WORD_MAX;
    break;
  case Affine_I:
    Atabcolumn[0].Rvalue = GT_WORD_MAX;
    Atabcolumn[0].Dvalue = GT_WORD_MAX;
    Atabcolumn[0].Ivalue = 0;
    break;
  default:
    Atabcolumn[0].Rvalue = 0;
    Atabcolumn[0].Dvalue = gap_opening;
    Atabcolumn[0].Ivalue = gap_opening;
  }
}

static void firstAtabRtabcolumn(GtAffinealignDPentry *Atabcolumn,
                                GtAffineAlignRtabentry *Rtabcolumn,
                                GtUword ulen,
                                GtUword gap_opening,
                                GtUword gap_extension,
                                GtAffineAlignEdge edge)
{
  GtUword rowindex;
  GtWord rdist, ddist,idist;
  firstAtabRtabentry(Atabcolumn, gap_opening, edge);

  Rtabcolumn[0].val_R.idx = 0;
  Rtabcolumn[0].val_D.idx = 0;
  Rtabcolumn[0].val_I.idx = 0;

  Rtabcolumn[0].val_R.edge = Affine_R;
  Rtabcolumn[0].val_D.edge = Affine_D;
  Rtabcolumn[0].val_I.edge = Affine_I;

  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    Atabcolumn[rowindex].Rvalue = GT_WORD_MAX;
    rdist = add_safe_max(Atabcolumn[rowindex-1].Rvalue,
                         gap_opening + gap_extension);
    ddist = add_safe_max(Atabcolumn[rowindex-1].Dvalue, gap_extension);
    idist = add_safe_max(Atabcolumn[rowindex-1].Dvalue,
                         gap_opening + gap_extension);
    Atabcolumn[rowindex].Dvalue = MIN3(rdist, ddist, idist);
    Atabcolumn[rowindex].Ivalue = GT_WORD_MAX;

    Atabcolumn[rowindex].Redge = Affine_X;
    Atabcolumn[rowindex].Dedge = gt_linearalign_affinegapcost_set_edge(rdist,
                                                                       ddist,
                                                                       idist);
    Atabcolumn[rowindex].Iedge = Affine_X;

    Rtabcolumn[rowindex].val_R.idx = rowindex;
    Rtabcolumn[rowindex].val_D.idx = rowindex;
    Rtabcolumn[rowindex].val_I.idx = rowindex;

    Rtabcolumn[rowindex].val_R.edge = Affine_R;
    Rtabcolumn[rowindex].val_D.edge = Affine_D;
    Rtabcolumn[rowindex].val_I.edge = Affine_I;
  }
}

static void nextAtabRtabcolumn(GtAffinealignDPentry *Atabcolumn,
                               GtAffineAlignRtabentry *Rtabcolumn,
                               const GtScoreHandler *scorehandler,
                               const GtUchar *useq,
                               GtUword ustart,
                               GtUword ulen,
                               GtUchar b,
                               GtUword midcolumn,
                               GtUword colindex)
{
  GtAffinealignDPentry northwestAffinealignDPentry, westAffinealignDPentry;
  GtAffineAlignRtabentry northwestRtabentry, westRtabentry;
  GtWord rowindex, rcost, rdist, ddist, idist, minvalue;
  GtUword gap_opening, gap_extension;

  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);

  northwestAffinealignDPentry = Atabcolumn[0];
  northwestRtabentry = Rtabcolumn[0];

  rdist = add_safe_max(Atabcolumn[0].Rvalue, gap_extension + gap_opening);
  ddist = add_safe_max(Atabcolumn[0].Dvalue, gap_extension + gap_opening);
  idist = add_safe_max(Atabcolumn[0].Ivalue, gap_extension);

  minvalue = MIN3(rdist, ddist, idist);
  Atabcolumn[0].Ivalue = minvalue;
  Atabcolumn[0].Rvalue = GT_WORD_MAX;
  Atabcolumn[0].Dvalue = GT_WORD_MAX;

  Atabcolumn[0].Redge = Affine_X;
  Atabcolumn[0].Dedge = Affine_X;
  Atabcolumn[0].Iedge = gt_linearalign_affinegapcost_set_edge(
                                                           rdist, ddist, idist);

  if (colindex > midcolumn)
  {
    northwestRtabentry = Rtabcolumn[0];
    Rtabcolumn[0].val_R.idx = Rtabcolumn[0].val_I.idx;
    Rtabcolumn[0].val_D.idx = Rtabcolumn[0].val_I.idx;
    Rtabcolumn[0].val_I.idx = Rtabcolumn[0].val_I.idx;

    Rtabcolumn[0].val_R.edge = Affine_X;
    Rtabcolumn[0].val_D.edge = Affine_X;
    Rtabcolumn[0].val_I.edge = Rtabcolumn[0].val_I.edge;
  }

  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    westAffinealignDPentry = Atabcolumn[rowindex];
    westRtabentry = Rtabcolumn[rowindex];

    rcost = gt_scorehandler_get_replacement(scorehandler,
                                            useq[ustart+rowindex-1], b);
    rdist = add_safe_max(northwestAffinealignDPentry.Rvalue, rcost);
    ddist = add_safe_max(northwestAffinealignDPentry.Dvalue, rcost);
    idist = add_safe_max(northwestAffinealignDPentry.Ivalue, rcost);

    minvalue = MIN3(rdist, ddist, idist);
    Atabcolumn[rowindex].Rvalue = minvalue;
    Atabcolumn[rowindex].Redge = gt_linearalign_affinegapcost_set_edge(rdist,
                                                                       ddist,
                                                                       idist);

    rdist = add_safe_max(Atabcolumn[rowindex-1].Rvalue,
                         gap_extension + gap_opening);
    ddist = add_safe_max(Atabcolumn[rowindex-1].Dvalue,gap_extension);
    idist = add_safe_max(Atabcolumn[rowindex-1].Ivalue,
                         gap_extension + gap_opening);

    minvalue = MIN3(rdist, ddist, idist);
    Atabcolumn[rowindex].Dvalue = minvalue;
    Atabcolumn[rowindex].Dedge = gt_linearalign_affinegapcost_set_edge(rdist,
                                                                       ddist,
                                                                       idist);

    rdist = add_safe_max(westAffinealignDPentry.Rvalue,
                         gap_extension + gap_opening);
    ddist = add_safe_max(westAffinealignDPentry.Dvalue,
                         gap_extension + gap_opening);
    idist = add_safe_max(westAffinealignDPentry.Ivalue, gap_extension);

    minvalue = MIN3(rdist, ddist, idist);
    Atabcolumn[rowindex].Ivalue = minvalue;
    Atabcolumn[rowindex].Iedge = gt_linearalign_affinegapcost_set_edge(rdist,
                                                                       ddist,
                                                                       idist);

    if (colindex > midcolumn)
    {
      Rtabcolumn[rowindex].val_R = get_Rtabentry(&northwestRtabentry,
                                                  Atabcolumn[rowindex].Redge);
      Rtabcolumn[rowindex].val_D = get_Rtabentry(&Rtabcolumn[rowindex-1],
                                                  Atabcolumn[rowindex].Dedge);
      Rtabcolumn[rowindex].val_I = get_Rtabentry(&westRtabentry,
                                                  Atabcolumn[rowindex].Iedge);
    }
    northwestAffinealignDPentry = westAffinealignDPentry;
    northwestRtabentry = westRtabentry;
  }
}

static GtUword evaluateallAtabRtabcolumns(GtAffinealignDPentry *Atabcolumn,
                                          GtAffineAlignRtabentry *Rtabcolumn,
                                          const GtScoreHandler *scorehandler,
                                          const GtUchar *useq,
                                          GtUword ustart,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vstart,
                                          GtUword vlen,
                                          GtUword midcolumn,
                                          GtAffineAlignEdge edge)
{
  GtUword colindex, gap_opening, gap_extension;

  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);

  firstAtabRtabcolumn(Atabcolumn, Rtabcolumn, ulen,
                      gap_opening, gap_extension, edge);

  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextAtabRtabcolumn(Atabcolumn,
                       Rtabcolumn,
                       scorehandler,
                       useq, ustart,ulen,
                       vseq[vstart+colindex-1],
                       midcolumn,
                       colindex);
  }

  return MIN3(Atabcolumn[ulen].Rvalue,
              Atabcolumn[ulen].Dvalue,
              Atabcolumn[ulen].Ivalue);
}

GtAffineAlignEdge gt_linearalign_affinegapcost_minAdditionalCosts(
                                              const GtAffinealignDPentry *entry,
                                              const GtAffineAlignEdge edge,
                                              GtUword gap_opening)
{
  GtUword rdist, ddist, idist;

  switch (edge) {
    case Affine_D:
      rdist = add_safe_max(entry->Rvalue, gap_opening);
      ddist = entry->Dvalue;
      idist = add_safe_max(entry->Ivalue, gap_opening);
     break;
    case Affine_I:
      rdist = add_safe_max(entry->Rvalue, gap_opening);
      ddist = add_safe_max(entry->Dvalue, gap_opening);
      idist = entry->Ivalue;
      break;
    default:
      rdist = entry->Rvalue;
      ddist = entry->Dvalue;
      idist = entry->Ivalue;
  }

  return gt_linearalign_affinegapcost_set_edge(rdist, ddist, idist);
}

#ifdef GT_THREADS_ENABLED
typedef struct{
  GtLinspaceManagement *spacemanager;
  const GtScoreHandler *scorehandler;
  const GtUchar *useq, * vseq;
  GtUword ustart, ulen, vstart, vlen,
          *Ctab, rowoffset, *threadcount;
  GtAffineAlignEdge from_edge, to_edge;
}GtAffineCrosspointthreadinfo;

static GtAffineCrosspointthreadinfo
              set_AffineCrosspointthreadinfo(GtLinspaceManagement *spacemanager,
                                             const GtScoreHandler *scorehandler,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen,
                                             GtUword *Ctab,
                                             GtUword rowoffset,
                                             GtAffineAlignEdge from_edge,
                                             GtAffineAlignEdge to_edge,
                                             GtUword *threadcount)
{
  GtAffineCrosspointthreadinfo threadinfo;
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
  threadinfo.from_edge = from_edge;
  threadinfo.to_edge = to_edge;
  threadinfo.threadcount = threadcount;

  return threadinfo;
}
static GtUword evaluateaffinecrosspoints(GtLinspaceManagement *spacemanager,
                                         const GtScoreHandler *scorehandler,
                                         const GtUchar *useq,
                                         GtUword ustart,
                                         GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart,
                                         GtUword vlen,
                                         GtUword *Ctab,
                                         GtUword rowoffset,
                                         GtAffineAlignEdge from_edge,
                                         GtAffineAlignEdge to_edge,
                                         GtUword *threadcount);

static void *evaluateaffinecrosspoints_thread_caller(void *data)
{
  GtAffineCrosspointthreadinfo *threadinfo =
                                         (GtAffineCrosspointthreadinfo *) data;
  (void) evaluateaffinecrosspoints(threadinfo->spacemanager,
                                   threadinfo->scorehandler,
                                   threadinfo->useq,
                                   threadinfo->ustart,
                                   threadinfo-> ulen,
                                   threadinfo->vseq,
                                   threadinfo->vstart,
                                   threadinfo->vlen,
                                   threadinfo->Ctab,
                                   threadinfo->rowoffset,
                                   threadinfo->from_edge,
                                   threadinfo->to_edge,
                                   threadinfo->threadcount);
  return NULL;
}
#endif

/* evaluate crosspoints in recursive way */
static GtUword evaluateaffinecrosspoints(GtLinspaceManagement *spacemanager,
                                         const GtScoreHandler *scorehandler,
                                         const GtUchar *useq,
                                         GtUword ustart,
                                         GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart,
                                         GtUword vlen,
                                         GtUword *Ctab,
                                         GtUword rowoffset,
                                         GtAffineAlignEdge from_edge,
                                         GtAffineAlignEdge to_edge,
                                         GT_UNUSED GtUword *threadcount)
{
  GtUword  midrow = 0, midcol = GT_DIV2(vlen), distance, colindex;
  GtAffineAlignEdge bottomtype, midtype = Affine_X;
  GtAffinealignDPentry *Atabcolumn = NULL;
  GtAffineAlignRtabentry *Rtabcolumn = NULL;

#ifdef GT_THREADS_ENABLED
  GtThread *t1 = NULL, *t2 = NULL;
  GtAffineCrosspointthreadinfo threadinfo1, threadinfo2;
#endif

  if (vlen >= 2UL)
  {
#ifdef GT_THREADS_ENABLED
    if (gt_jobs == 1)
    {
#endif
      if (gt_linspace_management_checksquare(spacemanager, ulen, vlen,
                                             sizeof (*Atabcolumn),
                                             sizeof (*Rtabcolumn)))
      {
        gt_affinealign_ctab(spacemanager, scorehandler, Ctab,
                            useq, ustart, ulen, vseq, vstart, vlen,
                            rowoffset, from_edge, to_edge);
        return 0;
      }
#ifdef GT_THREADS_ENABLED
   }
#endif
    Rtabcolumn = gt_linspace_management_get_rTabspace(spacemanager);
    Atabcolumn = gt_linspace_management_get_valueTabspace(spacemanager);
    Rtabcolumn = Rtabcolumn + rowoffset;
    Atabcolumn = Atabcolumn + rowoffset;

    distance = evaluateallAtabRtabcolumns(Atabcolumn,Rtabcolumn,
                                          scorehandler,
                                          useq, ustart, ulen,
                                          vseq, vstart, vlen,
                                          midcol, from_edge);

    bottomtype = gt_linearalign_affinegapcost_minAdditionalCosts(
                                 &Atabcolumn[ulen], to_edge,
                                 gt_scorehandler_get_gap_opening(scorehandler));
    switch (bottomtype) {
      case Affine_R:
        midrow = Rtabcolumn[ulen].val_R.idx;
        midtype = Rtabcolumn[ulen].val_R.edge;
        break;
      case Affine_D:
        midrow = Rtabcolumn[ulen].val_D.idx;
        midtype = Rtabcolumn[ulen].val_D.edge;
        break;
      case Affine_I:
        midrow = Rtabcolumn[ulen].val_I.idx;
        midtype = Rtabcolumn[ulen].val_I.edge;
        break;
      case Affine_X: /*never reach this line*/
        gt_assert(false);
    }
    Ctab[midcol] = rowoffset + midrow;
    gt_assert(midcol > 0);
    if (midrow == 0) {
      for (colindex = midcol-1; colindex > 0; colindex--)
        Ctab[colindex] = Ctab[midcol];
    }
    else{/* upper left corner */
      switch (midtype) {
        case Affine_R:
          if (midcol > 1)
            Ctab[midcol-1] = Ctab[midcol] == 0 ? 0: Ctab[midcol] - 1;

#ifdef GT_THREADS_ENABLED
          if (*threadcount + 1 > gt_jobs)
          {
#endif
            (void) evaluateaffinecrosspoints(spacemanager, scorehandler,
                                             useq, ustart, midrow-1,
                                             vseq, vstart, midcol-1,
                                             Ctab, rowoffset,
                                             from_edge, midtype,
                                             threadcount);
#ifdef GT_THREADS_ENABLED
          }
          else
          {
            threadinfo1=set_AffineCrosspointthreadinfo(spacemanager,
                                                       scorehandler,
                                                       useq, ustart, midrow-1,
                                                       vseq, vstart, midcol-1,
                                                       Ctab, rowoffset,
                                                       from_edge, midtype,
                                                       threadcount);
            (*threadcount)++;
            t1 = gt_thread_new(evaluateaffinecrosspoints_thread_caller,
                               &threadinfo1, NULL);
          }
#endif
          break;
        case Affine_D:
#ifdef GT_THREADS_ENABLED
          if (*threadcount + 1 > gt_jobs)
          {
#endif
          (void) evaluateaffinecrosspoints(spacemanager, scorehandler,
                                           useq, ustart, midrow-1,
                                           vseq, vstart, midcol,
                                           Ctab, rowoffset,
                                           from_edge,midtype,
                                           threadcount);
#ifdef GT_THREADS_ENABLED
         }
         else
         {
           threadinfo1=set_AffineCrosspointthreadinfo(spacemanager,
                                                      scorehandler,
                                                      useq, ustart, midrow-1,
                                                      vseq, vstart, midcol,
                                                      Ctab, rowoffset,
                                                      from_edge, midtype,
                                                      threadcount);

           (*threadcount)++;
           t1 = gt_thread_new(evaluateaffinecrosspoints_thread_caller,
                               &threadinfo1, NULL);
          }
#endif
          break;
        case Affine_I:
          if (midcol > 1)
            Ctab[midcol-1] = Ctab[midcol];
          (void) evaluateaffinecrosspoints(spacemanager, scorehandler,
                                           useq, ustart, midrow,
                                           vseq, vstart, midcol-1,
                                           Ctab, rowoffset,
                                           from_edge, midtype,
                                           threadcount);
          break;
        case Affine_X: /*never reach this line*/
                gt_assert(false);
      }
    }
   /*bottom right corner */
#ifdef GT_THREADS_ENABLED
    if (*threadcount + 1 > gt_jobs)
    {
#endif
      (void) evaluateaffinecrosspoints(spacemanager, scorehandler,
                                          useq, ustart+midrow, ulen-midrow,
                                          vseq, vstart+midcol, vlen-midcol,
                                          Ctab+midcol,rowoffset+midrow,
                                          midtype, to_edge, threadcount);
#ifdef GT_THREADS_ENABLED
    }
    else
    {
      threadinfo2 = set_AffineCrosspointthreadinfo(spacemanager, scorehandler,
                                                   useq, ustart + midrow,
                                                   ulen - midrow,
                                                   vseq, vstart + midcol,
                                                   vlen - midcol,
                                                   Ctab + midcol,
                                                   rowoffset + midrow,
                                                   midtype, to_edge,
                                                   threadcount);
      (*threadcount)++;
      t2 = gt_thread_new(evaluateaffinecrosspoints_thread_caller,
                         &threadinfo2, NULL);
    }

    if (t1 != NULL)
    {
      gt_thread_join(t1);
      (*threadcount)--;
      gt_thread_delete(t1);
    }
    if (t2 != NULL)
    {
      gt_thread_join(t2);
      (*threadcount)--;
      gt_thread_delete(t2);
    }

#endif
    return distance;
  }
  return 0;
}

static void affine_determineCtab0(GtUword *Ctab,
                                  GtLinspaceManagement *spacemanager,
                                  const GtScoreHandler *scorehandler,
                                  const GtUchar *useq,
                                  GtUword ustart,
                                  const GtUchar *vseq,
                                  GtUword vstart)
{
  GtAffinealignDPentry *Atabcolumn;

    if (Ctab[1]== 1 || Ctab[1] == 0)
      Ctab[0] = 0;
    else
    {
      gt_linspace_management_check(spacemanager,2*(Ctab[1]+1),Ctab[1],
                                   sizeof (*Atabcolumn),sizeof (Atabcolumn),0);
      /*gt_assert(vlen > 1);*/
      GtAffineAlignEdge to_edge_test = Affine_X;
      if (Ctab[1] == Ctab[2])
        to_edge_test = Affine_I;
      else
        to_edge_test = Affine_R;
      gt_affinealign_ctab(spacemanager, scorehandler, Ctab,
                          useq, ustart, Ctab[1], vseq, vstart,
                          1, 0, Affine_X, to_edge_test);
    }
}

/* calculating affine alignment in linear space */
GtUword gt_calc_affinealign_linear(GtLinspaceManagement *spacemanager,
                                   const GtScoreHandler *scorehandler,
                                   GtAlignment *align,
                                   const GtUchar *useq,
                                   GtUword ustart,
                                   GtUword ulen,
                                   const GtUchar *vseq,
                                   GtUword vstart,
                                   GtUword vlen)
{
  GtUword distance, *Ctab, threadcount = 1;
  GtWord gap_extension, gap_opening;
  GtAffinealignDPentry *Atabcolumn;
  GtAffineAlignRtabentry *Rtabcolumn;

  gt_linspace_management_set_ulen(spacemanager, ulen);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);
  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  if (ulen == 0UL)
  {
      distance = gt_reconstructalignment_trivial_insertion(align, vlen,
                                                           gap_extension);
      distance += gap_opening;
      return distance;
  }
  else if (vlen == 0UL)
  {
      distance = gt_reconstructalignment_trivial_deletion(align, ulen,
                                                          gap_extension);
      distance += gap_opening;
      return distance;
  }
  else if (vlen == 1UL)
  {
     gt_linspace_management_check(spacemanager, (ulen+1)*(vlen+1)-1, ulen,
                                  sizeof (*Atabcolumn), sizeof (Atabcolumn), 0);
    return gt_affinealign_with_Management(spacemanager, scorehandler, align,
                                   useq+ustart, ulen, vseq+vstart, vlen);
  }
  if (gt_linspace_management_checksquare(spacemanager, ulen, vlen,
                                         sizeof (*Atabcolumn),
                                         sizeof (*Rtabcolumn)))
  {
    return gt_affinealign_with_Management(spacemanager, scorehandler, align,
                                   useq+ustart, ulen, vseq+vstart, vlen);
  }
  else
  {
    gt_linspace_management_check(spacemanager, ulen, vlen, sizeof (*Atabcolumn),
                                 sizeof (*Rtabcolumn), sizeof (*Ctab));
    Ctab = gt_linspace_management_get_crosspointTabspace(spacemanager);
    Ctab[vlen] = ulen;
    distance = evaluateaffinecrosspoints(spacemanager, scorehandler,
                                         useq, ustart, ulen,
                                         vseq, vstart, vlen,
                                         Ctab, 0, Affine_X,
                                         Affine_X, &threadcount);

    affine_determineCtab0(Ctab, spacemanager, scorehandler,
                          useq, ustart, vseq, vstart);

    gt_reconstructalignment_from_Ctab(align, Ctab, useq, ustart, vseq,
                                      vstart, vlen, scorehandler);

  }
  return distance;
}

/* global alignment with affine gapcosts in linear space */
GtUword gt_linearalign_affinegapcost_compute_generic(
                                            GtLinspaceManagement *spacemanager,
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

  gt_alignment_set_seqs(align,useq+ustart, ulen, vseq+vstart, vlen);
  distance = gt_calc_affinealign_linear(spacemanager, scorehandler, align,
                                        useq, ustart, ulen,
                                        vseq, vstart, vlen);
  return distance;
}

/* global alignment with affine gapcosts in linear space
   with constant cost values*/
GtUword gt_linearalign_affinegapcost_compute(GtLinspaceManagement *spacemanager,
                                             GtAlignment *align,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen,
                                             GtUword matchcost,
                                             GtUword mismatchcost,
                                             GtUword gap_opening,
                                             GtUword gap_extension)
{
  GtUword distance;
  GtScoreHandler *scorehandler = gt_scorehandler_new(matchcost,
                                                     mismatchcost,
                                                     gap_opening,
                                                     gap_extension);

  gt_alignment_set_seqs(align,useq+ustart, ulen, vseq+vstart, vlen);
  distance = gt_linearalign_affinegapcost_compute_generic(spacemanager,
                                                          scorehandler, align,
                                                          useq, ustart, ulen,
                                                          vseq, vstart, vlen);
  gt_scorehandler_delete(scorehandler);
  return distance;
}

/*------------------------------local affine--------------------------------*/
static void firstAStabcolumn(GtAffinealignDPentry *Atabcolumn,
                             Starttabentry *Starttabcolumn,
                             const GtScoreHandler *scorehandler,
                             GtUword ulen)
{
  GtUword rowindex;
  GtWord gap_opening, gap_extension;

  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);

  Atabcolumn[0].Rvalue = GT_WORD_MIN;
  Atabcolumn[0].Dvalue = GT_WORD_MIN;
  Atabcolumn[0].Ivalue = GT_WORD_MIN;
  Atabcolumn[0].totalvalue = 0;

  Starttabcolumn[0].Rstart.a = 0;
  Starttabcolumn[0].Rstart.b = 0;
  Starttabcolumn[0].Dstart.a = 0;
  Starttabcolumn[0].Dstart.b = 0;
  Starttabcolumn[0].Istart.a = 0;
  Starttabcolumn[0].Istart.b = 0;

  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    Atabcolumn[rowindex].Rvalue = GT_WORD_MIN;
    Atabcolumn[rowindex].Dvalue = (gap_opening + gap_extension);
    Atabcolumn[rowindex].Ivalue = GT_WORD_MIN;
    Atabcolumn[rowindex].totalvalue = 0;

    Starttabcolumn[rowindex].Rstart.a = rowindex;
    Starttabcolumn[rowindex].Rstart.b = 0;
    Starttabcolumn[rowindex].Dstart.a = rowindex;
    Starttabcolumn[rowindex].Dstart.b = 0;
    Starttabcolumn[rowindex].Istart.a = rowindex;
    Starttabcolumn[rowindex].Istart.b = 0;
  }
}

static GtUwordPair setStarttabentry(GtWord entry, GtAffinealignDPentry *Atab,
                                    Starttabentry *Stab,
                                    GtWord replacement,
                                    GtWord gap_opening,
                                    GtWord gap_extension,
                                    const GtAffineAlignEdge edge)
{
  GtUwordPair start;
  switch (edge) {
    case Affine_R:
      if (entry == Atab->Rvalue + replacement)
         start = Stab->Rstart;
      else if (entry == Atab->Dvalue + replacement)
         start = Stab->Dstart;
      else if (entry == Atab->Ivalue + replacement)
         start = Stab->Istart;
      else
        start = Stab->Rstart;
      break;
    case Affine_D:
      if (entry == Atab->Rvalue + gap_opening + gap_extension)
         start = Stab->Rstart;
      else if (entry == Atab->Dvalue + gap_extension)
         start = Stab->Dstart;
      else if (entry == Atab->Ivalue + gap_opening + gap_extension)
         start = Stab->Istart;
      else
        start = Stab->Rstart;
      break;
    case Affine_I:
      if (entry == Atab->Rvalue + gap_opening + gap_extension)
         start = Stab->Rstart;
      else if (entry == Atab->Dvalue + gap_opening + gap_extension)
         start = Stab->Dstart;
      else if (entry == Atab->Ivalue + gap_extension)
         start = Stab->Istart;
      else
        start = Stab->Rstart;
      break;
    default:
      start.a = 0;
      start.b = 0;
  }
  return start;
}

static void nextAStabcolumn(GtAffinealignDPentry *Atabcolumn,
                            Starttabentry *Starttabcolumn,
                            const GtScoreHandler *scorehandler,
                            const GtUchar *useq, GtUword ustart,
                            GtUword ulen,
                            GtUchar b,
                            GtUword colindex,
                            GtMaxcoordvalue *max)
{
  GtAffinealignDPentry northwestAffinealignDPentry, westAffinealignDPentry;
  Starttabentry Snw, Swe;
  GtUword rowindex;
  GtWord gap_extension, gap_opening, replacement, temp, val1, val2;
  GtUwordPair start = {0};

  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);

  northwestAffinealignDPentry = Atabcolumn[0];
  Snw = Starttabcolumn[0];
  Atabcolumn[0].Rvalue = GT_WORD_MIN;
  Atabcolumn[0].Dvalue = GT_WORD_MIN;
  Atabcolumn[0].Ivalue = (gap_opening + gap_extension);
  temp = MAX3(Atabcolumn[0].Rvalue,
              Atabcolumn[0].Dvalue,
              Atabcolumn[0].Ivalue);
  Atabcolumn[0].totalvalue = ((temp > 0)? temp : 0);

  if (Atabcolumn[0].totalvalue == 0) {
    Starttabcolumn[0].Rstart.a = 0;
    Starttabcolumn[0].Rstart.b = colindex;
    Starttabcolumn[0].Dstart.a = 0;
    Starttabcolumn[0].Dstart.b = colindex;
    Starttabcolumn[0].Istart.a = 0;
    Starttabcolumn[0].Istart.b = colindex;
  }

  if (Atabcolumn[0].totalvalue > gt_maxcoordvalue_get_value(max))
    {
      if (Atabcolumn[0].totalvalue == Atabcolumn[0].Rvalue)
         start = Starttabcolumn[0].Rstart;
      else if (Atabcolumn[0].totalvalue == Atabcolumn[0].Dvalue)
         start = Starttabcolumn[0].Dstart;
      else if (Atabcolumn[0].totalvalue == Atabcolumn[0].Ivalue)
         start = Starttabcolumn[0].Istart;

      gt_maxcoordvalue_coord_update(max, Atabcolumn[0].totalvalue,
                                    start, 0, colindex);
    }
  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    westAffinealignDPentry = Atabcolumn[rowindex];
    Swe = Starttabcolumn[rowindex];

    /*calculate Rvalue*/
    replacement = gt_scorehandler_get_replacement(scorehandler,
                                                  useq[ustart+rowindex-1], b);

    Atabcolumn[rowindex].Rvalue =
              add_safe_min(northwestAffinealignDPentry.totalvalue, replacement);
    Starttabcolumn[rowindex].Rstart =
    setStarttabentry(Atabcolumn[rowindex].Rvalue, &northwestAffinealignDPentry,
                      &Snw, replacement,gap_opening, gap_extension, Affine_R);

    /*calculate Dvalue*/
    val1 = add_safe_min(Atabcolumn[rowindex-1].Dvalue, gap_extension);
    val2 = add_safe_min(Atabcolumn[rowindex-1].totalvalue,
                       (gap_opening+gap_extension));
    Atabcolumn[rowindex].Dvalue = MAX(val1,val2);
    Starttabcolumn[rowindex].Dstart =
    setStarttabentry(Atabcolumn[rowindex].Dvalue, &Atabcolumn[rowindex-1],
                     &Starttabcolumn[rowindex-1], replacement, gap_opening,
                     gap_extension,Affine_D);

    /*calculate Ivalue*/
    val1=(add_safe_min(westAffinealignDPentry.Ivalue,gap_extension));
    val2=(add_safe_min(westAffinealignDPentry.totalvalue,
                              gap_opening+gap_extension));
    Atabcolumn[rowindex].Ivalue = MAX(val1,val2);
    Starttabcolumn[rowindex].Istart =
    setStarttabentry(Atabcolumn[rowindex].Ivalue, &westAffinealignDPentry, &Swe,
                     replacement, gap_opening, gap_extension, Affine_I);

    /*calculate totalvalue*/
    temp = MAX3(Atabcolumn[rowindex].Rvalue,
                Atabcolumn[rowindex].Dvalue,
                Atabcolumn[rowindex].Ivalue);
    Atabcolumn[rowindex].totalvalue = temp > 0 ? temp : 0;

    /* set start indices for Atab-values*/
    if (Atabcolumn[rowindex].totalvalue == 0)
    {
      Starttabcolumn[rowindex].Rstart.a = rowindex;
      Starttabcolumn[rowindex].Rstart.b = colindex;
      Starttabcolumn[rowindex].Dstart.a = rowindex;
      Starttabcolumn[rowindex].Dstart.b = colindex;
      Starttabcolumn[rowindex].Istart.a = rowindex;
      Starttabcolumn[rowindex].Istart.b = colindex;
    }

    /*set new max*/
    if (Atabcolumn[rowindex].totalvalue > gt_maxcoordvalue_get_value(max))
    {
      if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Rvalue)
         start = Starttabcolumn[rowindex].Rstart;
      else if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Dvalue)
         start = Starttabcolumn[rowindex].Dstart;
      else if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Ivalue)
         start = Starttabcolumn[rowindex].Istart;

      gt_maxcoordvalue_coord_update(max, Atabcolumn[rowindex].totalvalue,
                                    start, rowindex, colindex);
    }
    northwestAffinealignDPentry = westAffinealignDPentry;
    Snw=Swe;
  }
}

static GtMaxcoordvalue *evaluateallAStabcolumns(GtLinspaceManagement *space,
                                            const GtScoreHandler *scorehandler,
                                            const GtUchar *useq,
                                            GtUword ustart,
                                            GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart,
                                            GtUword vlen)
{
  GtUword colindex;
  GtMaxcoordvalue *max;
  GtAffinealignDPentry *Atabcolumn;
  Starttabentry *Starttabcolumn;

  Atabcolumn = gt_linspace_management_get_valueTabspace(space);
  Starttabcolumn = gt_linspace_management_get_rTabspace(space);

  firstAStabcolumn(Atabcolumn, Starttabcolumn, scorehandler, ulen);

  max = gt_linspace_management_get_maxspace(space);
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextAStabcolumn(Atabcolumn, Starttabcolumn, scorehandler, useq, ustart,
                    ulen, vseq[vstart+colindex-1], colindex, max);
  }
  return max;
}

/* determining start and end of local alignment and call global function */
GtWord gt_linearalign_affinegapcost_compute_local_generic(
                                            GtLinspaceManagement *spacemanager,
                                            const GtScoreHandler *scorehandler,
                                            GtAlignment *align,
                                            const GtUchar *useq,
                                            GtUword ustart,
                                            GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart,
                                            GtUword vlen)
{
  GtUword ulen_part, ustart_part, vlen_part, vstart_part;
  GtWord score;
  GtAffinealignDPentry *Atabcolumn;
  Starttabentry *Starttabcolumn;
  GtMaxcoordvalue *max;

  gt_linspace_management_set_ulen(spacemanager, ulen);
  if (ulen == 0UL || vlen == 0UL)
  {
     /* empty alignment */
    return 0;
  }
  else if (vlen == 1UL)
  {
    gt_linspace_management_check_local(spacemanager,
                                       (ulen+1)*(vlen+1)-1, ulen,
                                       sizeof (*Atabcolumn),
                                       sizeof (Atabcolumn));
    return gt_affinealign_calculate_local_generic(spacemanager,scorehandler,
                                                  align,
                                                  useq, ustart, ulen,
                                                  vseq, vstart, vlen);
  }
  else if (gt_linspace_management_checksquare_local(spacemanager, ulen, vlen,
                                                    sizeof (*Atabcolumn),
                                                    sizeof (*Starttabcolumn)))
  {
    /* call alignment function for square space */
    return gt_affinealign_calculate_local_generic(spacemanager, scorehandler,
                                                  align, useq, ustart, ulen,
                                                  vseq, vstart, vlen);
  }

  gt_linspace_management_check_local(spacemanager, ulen, vlen,
                                     sizeof (*Atabcolumn),
                                     sizeof (*Starttabcolumn));

  max = evaluateallAStabcolumns(spacemanager, scorehandler,
                                useq, ustart, ulen,
                                vseq, vstart, vlen);

  score = gt_maxcoordvalue_get_value(max);

  if (gt_maxcoordvalue_get_length_safe(max))
  {
    GtScoreHandler *costhandler = gt_scorehandler2costhandler(scorehandler);

    ustart_part = ustart+(gt_maxcoordvalue_get_start(max)).a;
    vstart_part = vstart+(gt_maxcoordvalue_get_start(max)).b;
    ulen_part = gt_maxcoordvalue_get_row_length(max);
    vlen_part = gt_maxcoordvalue_get_col_length(max);

    gt_alignment_set_seqs(align,useq + ustart_part,ulen_part,
                                vseq + vstart_part,vlen_part);

    /* change score to cost, from maximizing to minimizing */

    gt_calc_affinealign_linear(spacemanager,
                               costhandler,
                               align,
                               useq, ustart_part, ulen_part,
                               vseq, vstart_part, vlen_part);
    gt_scorehandler_delete(costhandler);
  } else
  {
     /* empty alignment */
     return 0;
  }
  return score;
}

/* local alignment with linear gapcosts in linear space with constant costs*/
GtWord gt_linearalign_affinegapcost_compute_local(
                                             GtLinspaceManagement *spacemanager,
                                             GtAlignment *align,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen,
                                             GtWord matchscore,
                                             GtWord mismatchscore,
                                             GtWord gap_opening,
                                             GtWord gap_extension)
{
  GtWord score;

  GtScoreHandler *scorehandler = gt_scorehandler_new(matchscore,
                                                     mismatchscore,
                                                     gap_opening,
                                                     gap_extension);

  score = gt_linearalign_affinegapcost_compute_local_generic(spacemanager,
                                                            scorehandler, align,
                                                            useq, ustart, ulen,
                                                            vseq, vstart, vlen);
  gt_scorehandler_delete(scorehandler);
  return score;
}

/*----------------------------checkfunctions--------------------------*/
void gt_linearalign_affinegapcost_check(GT_UNUSED bool forward,
                                        const GtUchar *useq,
                                        GtUword ulen,
                                        const GtUchar *vseq,
                                        GtUword vlen)
{
  GtAlignment *align;
  GtUword affine_score1, affine_score2, affine_score3,
          matchcost = 0, mismatchcost = 4, gap_opening = 4, gap_extension = 1;
  GtLinspaceManagement *spacemanager;
  GtScoreHandler *scorehandler;

  gt_assert(useq && vseq);
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

  /* affinealign (square) handles lower/upper cases  in another way*/
  scorehandler = gt_scorehandler_new(matchcost, mismatchcost,
                                     gap_opening, gap_extension);
  gt_scorehandler_plain(scorehandler);
  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);
  spacemanager = gt_linspace_management_new();

  affine_score1 = gt_calc_affinealign_linear(spacemanager, scorehandler, align,
                                             useq, 0, ulen,
                                             vseq, 0, vlen);
  gt_linspace_management_delete(spacemanager);
  affine_score2 = gt_alignment_eval_with_affine_score(align,  true,matchcost,
                                                      mismatchcost, gap_opening,
                                                      gap_extension);
  gt_alignment_delete(align);
  gt_scorehandler_delete(scorehandler);
  if (affine_score1 != affine_score2)
  {
    fprintf(stderr,"gt_calc_affinealign_linear = "GT_WU" != "GT_WU
            " = gt_alignment_eval_with_affine_score\n", affine_score1,
                                                        affine_score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  align = gt_affinealign(useq, ulen, vseq, vlen, matchcost,
                         mismatchcost, gap_opening, gap_extension);

  affine_score3 = gt_alignment_eval_with_affine_score(align,  true, matchcost,
                                                      mismatchcost, gap_opening,
                                                      gap_extension);

  if (affine_score1 != affine_score3)
  {
    fprintf(stderr,"gt_calc_affinealign_linear = "GT_WU" != "GT_WU
            " = gt_affinealign\n", affine_score1, affine_score3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_alignment_delete(align);
}

void gt_linearalign_affinegapcost_check_local(GT_UNUSED bool forward,
                                              const GtUchar *useq,
                                              GtUword ulen,
                                              const GtUchar *vseq,
                                              GtUword vlen)
{
  GtAlignment *align;
  GtWord affine_score1, affine_score2, affine_score3, affine_score4,
         matchscore = 6, mismatchscore = -3,
         gap_opening = -2, gap_extension = -1;
  GtScoreHandler *scorehandler;
  GtLinspaceManagement *spacemanager;

  gt_assert(useq && vseq);
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

  scorehandler = gt_scorehandler_new(matchscore, mismatchscore,
                                     gap_opening, gap_extension);
  align = gt_alignment_new();
  spacemanager = gt_linspace_management_new();
  affine_score1 = gt_linearalign_affinegapcost_compute_local_generic(
                                                            spacemanager,
                                                            scorehandler, align,
                                                            useq, 0, ulen,
                                                            vseq,  0, vlen);
  gt_linspace_management_delete(spacemanager);
  gt_scorehandler_delete(scorehandler);

  affine_score2 = gt_alignment_eval_with_affine_score(align, true, matchscore,
                                    mismatchscore, gap_opening, gap_extension);

  if (affine_score1 != affine_score2)
  {
    fprintf(stderr,"gt_linearalign_affinegapcost_compute_local_generic ="
            " "GT_WD"!= "GT_WD
            " = gt_alignment_eval_with_affine_score\n", affine_score1,
                                                        affine_score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_alignment_reset(align);
  affine_score3 = gt_affinealign_calculate_local(NULL, align, useq, 0,
                                                 ulen, vseq, 0, vlen,
                                                 matchscore, mismatchscore,
                                                 gap_opening, gap_extension);

  if (affine_score1 != affine_score3)
  {
    fprintf(stderr,"gt_calc_affinealign_linear_local = "GT_WD" != "GT_WD
            " = affinealign_in_square_space_local\n", affine_score1,
                                                        affine_score3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  affine_score4 = gt_alignment_eval_with_affine_score(align, true, matchscore,
                                    mismatchscore, gap_opening, gap_extension);

  if (affine_score3 != affine_score4)
  {
    fprintf(stderr,"affinealign_in_square_space_local = "GT_WD" != "GT_WD
            " = gt_alignment_eval_generic_with_affine_score\n", affine_score3,
                                                        affine_score4);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_alignment_delete(align);
}
