/*
  Copyright (C) 2015 Annika Seidel, <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2007-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2015 Center for Bioinformatics, University of Hamburg

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
#include <limits.h>
#include "core/assert_api.h"
#include "core/array2dim_api.h"
#include "core/minmax.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/maxcoordvalue.h"

#include "extended/affinealign.h"

static void affinealign_fill_table(GtAffinealignDPentry **dptable,
                                   const GtUchar *u, GtUword ulen,
                                   const GtUchar *v, GtUword vlen,
                                   GtUword matchcost, GtUword mismatchcost,
                                   GtUword gap_opening, GtUword gap_extension,
                                   GtAffineAlignEdge edge,
                                   const GtScoreHandler *scorehandler)
{
  GtUword i, j, Rvalue, Dvalue, Ivalue, minvalue,rcost;
  gt_assert(dptable && u && v);

  if (scorehandler != NULL)/*else work with given constant values*/
  {
    gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
    gap_extension = gt_scorehandler_get_gapscore(scorehandler);
  }

  for (i = 0; i <= ulen; i++) {
    for (j = 0; j <= vlen; j++) {
      if (!i && !j) {
        switch (edge) {
          case Affine_R:
            dptable[0][0].Rvalue = 0;
            dptable[0][0].Dvalue = GT_WORD_MAX;
            dptable[0][0].Ivalue = GT_WORD_MAX;
            break;
          case Affine_D:
            dptable[0][0].Rvalue = GT_WORD_MAX;
            dptable[0][0].Dvalue = 0;
            dptable[0][0].Ivalue = GT_WORD_MAX;
            break;
          case Affine_I:
            dptable[0][0].Rvalue = GT_WORD_MAX;
            dptable[0][0].Dvalue = GT_WORD_MAX;
            dptable[0][0].Ivalue = 0;
            break;
          default:
            dptable[0][0].Rvalue = 0;
            dptable[0][0].Dvalue = gap_opening;
            dptable[0][0].Ivalue = gap_opening;
        }
      }
      else {
        /* compute A_affine(i,j,R) */
        if (!i || !j)
          dptable[i][j].Rvalue = LONG_MAX;
        else {
          if (scorehandler != NULL)
            rcost = gt_scorehandler_get_replacement(scorehandler,
                                                    u[i-1], v[j-1]);
          else
            rcost  = (u[i-1] == v[j-1]) ? matchcost : mismatchcost;

          Rvalue = add_safe_max(dptable[i-1][j-1].Rvalue, rcost);
          Dvalue = add_safe_max(dptable[i-1][j-1].Dvalue, rcost);
          Ivalue = add_safe_max(dptable[i-1][j-1].Ivalue, rcost);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          gt_assert(minvalue != LONG_MAX);
          dptable[i][j].Rvalue = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Redge = Affine_R;
          else if (Dvalue == minvalue)
            dptable[i][j].Redge = Affine_D;
          else /* Ivalue == minvalue */
            dptable[i][j].Redge = Affine_I;
        }
        /* compute A_affine(i,j,D) */
        if (!i)
          dptable[i][j].Dvalue = LONG_MAX;
        else {
          Rvalue = add_safe_max(dptable[i-1][j].Rvalue,
                                gap_opening + gap_extension);
          Dvalue = add_safe_max(dptable[i-1][j].Dvalue, gap_extension);
          Ivalue = add_safe_max(dptable[i-1][j].Ivalue,
                                gap_opening + gap_extension);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          gt_assert(minvalue != ULONG_MAX);
          dptable[i][j].Dvalue = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Dedge = Affine_R;
          else if (Dvalue == minvalue)
            dptable[i][j].Dedge = Affine_D;
          else /* Ivalue == minvalue */
            dptable[i][j].Dedge = Affine_I;
        }
        /* compute A_affine(i,j,I) */
        if (!j)
          dptable[i][j].Ivalue = LONG_MAX;
        else {
          Rvalue = add_safe_max(dptable[i][j-1].Rvalue,
                                gap_opening + gap_extension);
          Dvalue = add_safe_max(dptable[i][j-1].Dvalue,
                                gap_opening + gap_extension);
          Ivalue = add_safe_max(dptable[i][j-1].Ivalue, gap_extension);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          gt_assert(minvalue != LONG_MAX);
          dptable[i][j].Ivalue = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Iedge = Affine_R;
          else if (Dvalue == minvalue)
            dptable[i][j].Iedge = Affine_D;
          else /* Ivalue == minvalue */
            dptable[i][j].Iedge = Affine_I;
        }
      }
    }
  }
}

GtWord gt_affinealign_traceback(GtAlignment *align,
                                GtAffinealignDPentry * const *dptable,
                                GtUword i, GtUword j)
{
  GtWord minvalue;
  GtAffineAlignEdge edge;
  gt_assert(align && dptable);
  /* determine min{A_affine(m,n,x) | x in {R,D,I}} */
  minvalue = MIN3(dptable[i][j].Rvalue, dptable[i][j].Dvalue,
                  dptable[i][j].Ivalue);
  if (dptable[i][j].Rvalue == minvalue)
    edge = Affine_R;
  else if (dptable[i][j].Dvalue == minvalue)
    edge = Affine_D;
  else /* dptable[i][j].Ivalue == minvalue */
    edge = Affine_I;
  /* backtracing */
  while (i > 0 || j > 0) {
    switch (edge) {
      case Affine_R:
        gt_assert(dptable[i][j].Rvalue != LONG_MAX);
        gt_alignment_add_replacement(align);
        edge = dptable[i][j].Redge;
        gt_assert(i > 0 && j > 0);
        i--;
        j--;
        break;
      case Affine_D:
        gt_alignment_add_deletion(align);
        edge = dptable[i][j].Dedge;
        gt_assert(i);
        i--;
        break;
      case Affine_I:
        gt_alignment_add_insertion(align);
        edge = dptable[i][j].Iedge;
        gt_assert(j);
        j--;
        break;
      default:
        gt_assert(false);
    }
  }
  return minvalue;
}

GtAlignment* gt_affinealign(const GtUchar *u, GtUword ulen,
                            const GtUchar *v, GtUword vlen,
                            GtUword matchcost, GtUword mismatchcost,
                            GtUword gap_opening,
                            GtUword gap_extension)
{
  GtAffinealignDPentry **dptable;
  GtAlignment *align;

  gt_assert(u && v);
  gt_array2dim_malloc(dptable, ulen+1, vlen+1);
  affinealign_fill_table(dptable, u, ulen, v, vlen, matchcost, mismatchcost,
                         gap_opening, gap_extension, Affine_X, NULL);
  align = gt_alignment_new_with_seqs(u, ulen,  v, vlen);
  gt_affinealign_traceback(align, dptable, ulen, vlen);
  gt_array2dim_delete(dptable);
  return align;
}

GtWord gt_affinealign_with_Management(GtLinspaceManagement *spacemanager,
                                     const GtScoreHandler *scorehandler,
                                     GtAlignment *align,
                                     const GtUchar *u, GtUword ulen,
                                     const GtUchar *v, GtUword vlen)
{
  GtAffinealignDPentry **dptable;
  GtUword idx;
  gt_assert(u && v && spacemanager && scorehandler);

  gt_assert((ulen+1)*(vlen+1)*sizeof(**dptable) <=
             gt_linspace_management_get_valueTabsize(spacemanager));

  dptable = gt_linspace_management_get_rTabspace(spacemanager);
  *dptable = gt_linspace_management_get_valueTabspace(spacemanager);

  for (idx = 1; idx < ulen + 1; idx++)
    dptable[idx] = dptable[idx-1] + vlen + 1;

  affinealign_fill_table(dptable, u, ulen, v, vlen, GT_WORD_MAX, GT_WORD_MAX,
                         GT_WORD_MAX, GT_WORD_MAX, Affine_X, scorehandler);
  return gt_affinealign_traceback(align, dptable, ulen, vlen);
}

static void evaluate_affinecrosspoints_from_2dimtab(GtUword *Ctab,
                                            GtAffinealignDPentry **Atabcolumn,
                                            const GtScoreHandler *scorehandler,
                                            GtUword ulen, GtUword vlen,
                                            GtUword rowoffset,
                                            GtAffineAlignEdge edge)
{
  GtUword i, j, gap_opening;
  gt_assert(Atabcolumn != NULL);
  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);

  i = ulen;
  j = vlen;
  edge = gt_linearalign_affinegapcost_minAdditionalCosts(&Atabcolumn[i][j],
                                                         edge, gap_opening);

  while (i > 0 || j > 1) {
    switch (edge) {
      case Affine_R:
        gt_assert(Atabcolumn[i][j].Rvalue != LONG_MAX);
        Ctab[j-1] = i-1 + rowoffset;
        edge = Atabcolumn[i][j].Redge;
        gt_assert(i > 0 && j > 0);
        i--;
        j--;
        break;
      case Affine_D:
        edge = Atabcolumn[i][j].Dedge;
        gt_assert(i);
        i--;
        break;
      case Affine_I:
        Ctab[j-1] = i + rowoffset;
        edge = Atabcolumn[i][j].Iedge;
        gt_assert(j);
        j--;
        break;
      default:
        gt_assert(false);
    }
  }
}

void gt_affinealign_ctab(GtLinspaceManagement *spacemanager,
                         const GtScoreHandler *scorehandler,
                         GtUword *Ctab,
                         const GtUchar *useq,
                         GtUword ustart,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vstart,
                         GtUword vlen,
                         GtUword rowoffset,
                         GtAffineAlignEdge from_edge,
                         GtAffineAlignEdge to_edge)
{
  GtAffinealignDPentry **dptable;
  GtUword idx;
  gt_assert(Ctab && spacemanager && scorehandler);
  gt_assert((ulen+1)*(vlen+1)*sizeof(**dptable) <=
             gt_linspace_management_get_valueTabsize(spacemanager));

  dptable = gt_linspace_management_get_rTabspace(spacemanager);
  *dptable = gt_linspace_management_get_valueTabspace(spacemanager);

  for (idx = 1; idx < ulen+1; idx++)
    dptable[idx] = dptable[idx-1] + vlen + 1;

  affinealign_fill_table(dptable, &useq[ustart], ulen, &vseq[vstart], vlen,
                         GT_WORD_MAX, GT_WORD_MAX, GT_WORD_MAX,
                         GT_WORD_MAX, from_edge, scorehandler);

  evaluate_affinecrosspoints_from_2dimtab(Ctab, dptable, scorehandler, ulen,
                                          vlen, rowoffset, to_edge);

}

/* local */
static GtWord affinealign_fill_table_local(GtAffinealignDPentry **Atabcolumn,
                                           const GtScoreHandler *scorehandler,
                                           GtMaxcoordvalue *max,
                                           const GtUchar *useq, GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq, GtUword vstart,
                                           GtUword vlen,
                                           GtWord matchscore,
                                           GtWord mismatchscore,
                                           GtWord gap_opening,
                                           GtWord gap_extension)
{
  GtUword i,j;
  GtWord Rvalue, Dvalue, Ivalue, totalvalue, replacement, temp;

  if (scorehandler != NULL)/*else work with given constant values*/
  {
    gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
    gap_extension = gt_scorehandler_get_gapscore(scorehandler);
  }
  Atabcolumn[0][0].Rvalue = GT_WORD_MIN;
  Atabcolumn[0][0].Dvalue = GT_WORD_MIN;
  Atabcolumn[0][0].Ivalue = GT_WORD_MIN;
  Atabcolumn[0][0].totalvalue = 0;
  /* set backtracing edge */
  Atabcolumn[0][0].Redge = Affine_X;
  Atabcolumn[0][0].Dedge = Affine_X;
  Atabcolumn[0][0].Iedge = Affine_X;

  for (i = 1; i <= ulen; i++)
  {
    Atabcolumn[i][0].Rvalue = GT_WORD_MIN;
    Atabcolumn[i][0].Dvalue = gap_opening + gap_extension;
    Atabcolumn[i][0].Ivalue = GT_WORD_MIN;
    Atabcolumn[i][0].totalvalue = 0;
    /* set backtracing edge */
    Atabcolumn[i][0].Redge = Affine_X;
    Atabcolumn[i][0].Dedge = Affine_D;
    Atabcolumn[i][0].Iedge = Affine_X;
  }

  for (j = 1UL; j <= vlen; j++)
  {
    Atabcolumn[0][j].Rvalue = GT_WORD_MIN;
    Atabcolumn[0][j].Dvalue = GT_WORD_MIN;
    Atabcolumn[0][j].Ivalue = gap_opening + gap_extension;
    temp = MAX3(Atabcolumn[0][j].Rvalue,
                Atabcolumn[0][j].Dvalue,
                Atabcolumn[0][j].Ivalue);
    Atabcolumn[0][j].totalvalue = ((temp > 0)? temp : 0);
    /* set backtracing edge */
    Atabcolumn[0][j].Redge = Affine_X;
    Atabcolumn[0][j].Dedge = Affine_X;
    Atabcolumn[0][j].Iedge = Affine_I;

    if (Atabcolumn[0][j].totalvalue > gt_maxcoordvalue_get_value(max))
    {
      gt_maxcoordvalue_coord_update_without_start(max,
                                                  Atabcolumn[0][j].totalvalue,
                                                  0, j);
    }
    for (i = 1; i <= ulen; i++)
    {

      /*calculate Rvalue*/
      if (scorehandler != NULL)
      {
        replacement = gt_scorehandler_get_replacement(scorehandler,
                                useq[ustart+i-1], vseq[vstart+j-1]);
      }
      else
      {
        replacement = (tolower((int)useq[ustart+i-1]) ==
                       tolower((int)vseq[vstart+j-1]) ?
                       matchscore : mismatchscore);
      }
      Rvalue = add_safe_min(Atabcolumn[i-1][j-1].Rvalue, replacement);
      Dvalue = add_safe_min(Atabcolumn[i-1][j-1].Dvalue, replacement);
      Ivalue = add_safe_min(Atabcolumn[i-1][j-1].Ivalue, replacement);
      totalvalue = add_safe_min(Atabcolumn[i-1][j-1].totalvalue, replacement);
      Atabcolumn[i][j].Rvalue = MAX(MAX(Rvalue, Dvalue),MAX(Ivalue,totalvalue));
      /* set backtracing edge */
      if (Rvalue == Atabcolumn[i][j].Rvalue)
        Atabcolumn[i][j].Redge = Affine_R;
      else if (Dvalue == Atabcolumn[i][j].Rvalue)
        Atabcolumn[i][j].Redge = Affine_D;
      else if (Ivalue == Atabcolumn[i][j].Rvalue)
        Atabcolumn[i][j].Redge = Affine_I;
      else
        Atabcolumn[i][j].Redge = Affine_X;

      /*calculate Dvalue*/
      Rvalue = add_safe_min(Atabcolumn[i-1][j].Rvalue,
                            gap_opening + gap_extension);
      Dvalue = add_safe_min(Atabcolumn[i-1][j].Dvalue, gap_extension);
      Ivalue = add_safe_min(Atabcolumn[i-1][j].Ivalue,
                            gap_opening + gap_extension);
      totalvalue = add_safe_min(Atabcolumn[i-1][j].totalvalue,
                                gap_opening + gap_extension);
      Atabcolumn[i][j].Dvalue = MAX(MAX(Rvalue, Dvalue),MAX(Ivalue,totalvalue));
      /* set backtracing edge */
      if (Rvalue == Atabcolumn[i][j].Dvalue)
        Atabcolumn[i][j].Dedge = Affine_R;
      else if (Dvalue == Atabcolumn[i][j].Dvalue)
        Atabcolumn[i][j].Dedge = Affine_D;
      else if (Ivalue == Atabcolumn[i][j].Dvalue)
        Atabcolumn[i][j].Dedge = Affine_I;
      else
        Atabcolumn[i][j].Dedge = Affine_X;

      /*calculate Ivalue*/
      Rvalue = add_safe_min(Atabcolumn[i][j-1].Rvalue,
                            gap_extension + gap_opening);
      Dvalue = add_safe_min(Atabcolumn[i][j-1].Dvalue,
                            gap_extension + gap_opening);
      Ivalue = add_safe_min(Atabcolumn[i][j-1].Ivalue,gap_extension);
      totalvalue = add_safe_min(Atabcolumn[i][j-1].totalvalue,
                                gap_extension + gap_opening);
      Atabcolumn[i][j].Ivalue = MAX(MAX(Rvalue, Dvalue),MAX(Ivalue,totalvalue));
      /* set backtracing edge */
      if (Rvalue == Atabcolumn[i][j].Ivalue)
        Atabcolumn[i][j].Iedge = Affine_R;
      else if (Dvalue == Atabcolumn[i][j].Ivalue)
        Atabcolumn[i][j].Iedge = Affine_D;
      else if (Ivalue == Atabcolumn[i][j].Ivalue)
        Atabcolumn[i][j].Iedge = Affine_I;
      else
        Atabcolumn[i][j].Iedge = Affine_X;

      /*calculate totalvalue*/
      temp = MAX3(Atabcolumn[i][j].Rvalue,
                  Atabcolumn[i][j].Dvalue,
                  Atabcolumn[i][j].Ivalue);
      Atabcolumn[i][j].totalvalue = temp > 0 ? temp : 0;

      /*set new max*/
      if (Atabcolumn[i][j].totalvalue > gt_maxcoordvalue_get_value(max))
      {
        gt_maxcoordvalue_coord_update_without_start(max,
                                                    Atabcolumn[i][j].totalvalue,
                                                    i,j);
      }
    }
  }

  return gt_maxcoordvalue_get_value(max);
}

static void affinealign_traceback_local(GtAlignment *align,
                                        GtAffinealignDPentry **dptable,
                                        GtMaxcoordvalue *max)
{
  GtWord maxvalue;
  GtUword i,j;
  GtUwordPair max_end;
  GtAffineAlignEdge edge;
  gt_assert(align && dptable);

  max_end = gt_maxcoordvalue_get_end(max);
  i = max_end.a;
  j = max_end.b;

  maxvalue = MAX(MAX(dptable[i][j].Rvalue, dptable[i][j].Dvalue),
                  MAX(dptable[i][j].Ivalue,dptable[i][j].totalvalue));
  gt_assert(gt_maxcoordvalue_get_value(max) == maxvalue);

  if (dptable[i][j].Rvalue == maxvalue)
    edge = Affine_R;
  else if (dptable[i][j].Dvalue == maxvalue)
    edge = Affine_D;
  else /* dptable[i][j].Ivalue == maxvalue */
    edge = Affine_I;

  /* backtracing */
  while (edge != Affine_X && (i > 0 || j > 0)) {
    switch (edge) {
      case Affine_R:
        gt_assert(dptable[i][j].Rvalue != LONG_MAX);
        gt_alignment_add_replacement(align);
        edge = dptable[i][j].Redge;
        gt_assert(i > 0 && j > 0);
        i--;
        j--;
        break;
      case Affine_D:
        gt_alignment_add_deletion(align);
        edge = dptable[i][j].Dedge;
        gt_assert(i);
        i--;
        break;
      case Affine_I:
        gt_alignment_add_insertion(align);
        edge = dptable[i][j].Iedge;
        gt_assert(j);
        j--;
        break;
      default:
        gt_assert(false);
    }
  }
  gt_maxcoordvalue_set_start(max,i,j);
}

GtWord gt_affinealign_calculate_local_generic(GtLinspaceManagement
                                             *spacemanager,
                                             const GtScoreHandler *scorehandler,
                                             GtAlignment *align,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen)
{
  GtWord score = 0;
  GtAffinealignDPentry **Atabcolumn;
  GtMaxcoordvalue *max;

  gt_assert(align != NULL);
  if (spacemanager == NULL)
  {
    /*use it in normally case*/
    gt_array2dim_malloc(Atabcolumn,(ulen+1),(vlen+1));
    max = gt_maxcoordvalue_new();
  }
  else
  {
    /*use it in lineraspace context*/
    gt_assert((ulen+1)*(vlen+1)*sizeof(**Atabcolumn) <=
               gt_linspace_management_get_valueTabsize(spacemanager));

    Atabcolumn = gt_linspace_management_get_rTabspace(spacemanager);
    *Atabcolumn = gt_linspace_management_get_valueTabspace(spacemanager);

    GtUword idx;
    for (idx = 1;idx < ulen+1;idx++)
      Atabcolumn[idx] = Atabcolumn[idx-1] + vlen + 1;
    max = gt_linspace_management_get_maxspace(spacemanager);
  }

  score = affinealign_fill_table_local(Atabcolumn, scorehandler, max,
                                       useq, ustart, ulen,
                                       vseq, vstart, vlen,
                                       GT_WORD_MAX, GT_WORD_MAX,
                                       GT_WORD_MAX, GT_WORD_MAX);

  /* reconstruct local alignment from 2dimarray Atabcolumn */
  affinealign_traceback_local(align, Atabcolumn, max);

  if (gt_maxcoordvalue_get_length_safe(max))
  {
    ustart = ustart + (gt_maxcoordvalue_get_start(max)).a;
    vstart = vstart + (gt_maxcoordvalue_get_start(max)).b;
    ulen = gt_maxcoordvalue_get_row_length(max);
    vlen = gt_maxcoordvalue_get_col_length(max);

    gt_alignment_set_seqs(align, &useq[ustart], ulen,
                                 &vseq[vstart], vlen);
  }

  if (spacemanager == NULL)
  {
    gt_array2dim_delete(Atabcolumn);
    gt_maxcoordvalue_delete(max);
  }
  return score;
}

GtWord gt_affinealign_calculate_local(GtLinspaceManagement *spacemanager,
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
  score = gt_affinealign_calculate_local_generic(spacemanager, scorehandler,
                                                 align, useq, ustart, ulen,
                                                 vseq, vstart, vlen);
  gt_scorehandler_delete(scorehandler);
  return score;
}
