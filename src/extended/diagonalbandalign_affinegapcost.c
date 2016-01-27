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
#include "core/array2dim_api.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/divmodmul.h"
#include "core/ma_api.h"
#include "extended/affinealign.h"
#include "extended/diagonalbandalign.h"
#include "extended/diagonalbandalign_affinegapcost.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/linspace_management.h"
#include "extended/reconstructalignment.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

static void diagonalband_fillDPtab_affine(GtAffinealignDPentry **Atabcolumn,
                                          const GtUchar *useq,
                                          GtUword ustart,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vstart,
                                          GtUword vlen,
                                          GtWord left_dist,
                                          GtWord right_dist,
                                          GtAffineAlignEdge from_edge,
                                          GtAffineAlignEdge edge,
                                          const GtScoreHandler *scorehandler)
{
  GtUword i,j, low_row, high_row, gap_opening, gap_extension;
  GtWord rcost, r_dist, d_dist, i_dist, minvalue;

  gt_assert(Atabcolumn && scorehandler);
  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }

  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);
  low_row = 0;
  high_row = -left_dist;

  /* first entry */
   switch (edge) {
    case Affine_R:
      Atabcolumn[0][0].Rvalue = 0;
      Atabcolumn[0][0].Redge = from_edge;
      Atabcolumn[0][0].Dvalue = GT_WORD_MAX;
      Atabcolumn[0][0].Ivalue = GT_WORD_MAX;
      break;
    case Affine_D:
      Atabcolumn[0][0].Rvalue = GT_WORD_MAX;
      Atabcolumn[0][0].Dvalue = 0;
      Atabcolumn[0][0].Dedge = from_edge;
      Atabcolumn[0][0].Ivalue = GT_WORD_MAX;
      break;
    case Affine_I:
      Atabcolumn[0][0].Rvalue = GT_WORD_MAX;
      Atabcolumn[0][0].Dvalue = GT_WORD_MAX;
      Atabcolumn[0][0].Ivalue = 0;
      Atabcolumn[0][0].Iedge = from_edge;
      break;
    default:
      Atabcolumn[0][0].Rvalue = 0;
      Atabcolumn[0][0].Dvalue = gap_opening;
      Atabcolumn[0][0].Ivalue = gap_opening;
    }
  /* first column */
  for (i = 1; i <= high_row; i++)
  {
    Atabcolumn[i][0].Rvalue = GT_WORD_MAX;
    r_dist = add_safe_max(Atabcolumn[i-1][0].Rvalue,
                         gap_opening + gap_extension);
    d_dist = add_safe_max(Atabcolumn[i-1][0].Dvalue, gap_extension);
    i_dist = add_safe_max(Atabcolumn[i-1][0].Ivalue,
                         gap_opening + gap_extension);
    Atabcolumn[i][0].Dvalue = MIN3(r_dist, d_dist, i_dist);
    Atabcolumn[i][0].Ivalue = GT_WORD_MAX;

    Atabcolumn[i][0].Redge = Affine_X;
    Atabcolumn[i][0].Dedge = gt_linearalign_affinegapcost_set_edge(r_dist,
                                                                   d_dist,
                                                                   i_dist);
    Atabcolumn[i][0].Iedge = Affine_X;
  }
  for (; i <= ulen; i++)
  {
    /* invalid values */
    Atabcolumn[i][0].Rvalue = GT_WORD_MAX;
    Atabcolumn[i][0].Dvalue = GT_WORD_MAX;
    Atabcolumn[i][0].Ivalue = GT_WORD_MAX;
  }

  /* next columns */
  for (j = 1; j <= vlen; j++)
  {
    /* below diagonal band*/
    for (i = 0; i <= low_row; i++)
    {
      if (j <= right_dist)
      {
        Atabcolumn[i][j].Redge = Affine_X;
        Atabcolumn[i][j].Dedge = Affine_X;
        r_dist = add_safe_max(Atabcolumn[i][j-1].Rvalue,
                               gap_extension + gap_opening);
        d_dist = add_safe_max(Atabcolumn[i][j-1].Dvalue,
                               gap_extension + gap_opening);
        i_dist = add_safe_max(Atabcolumn[i][j-1].Ivalue,gap_extension);

        minvalue = MIN3(r_dist, d_dist, i_dist);
        Atabcolumn[i][j].Ivalue = minvalue;
        Atabcolumn[i][j].Rvalue = GT_WORD_MAX;
        Atabcolumn[i][j].Dvalue = GT_WORD_MAX;

        Atabcolumn[i][j].Iedge = gt_linearalign_affinegapcost_set_edge(r_dist,
                                                         d_dist,
                                                         i_dist);
      }
      else{
        Atabcolumn[i][j].Rvalue = GT_WORD_MAX;
        Atabcolumn[i][j].Dvalue = GT_WORD_MAX;
        Atabcolumn[i][j].Ivalue = GT_WORD_MAX;
        Atabcolumn[i][j].Iedge = Affine_X;
      }
    }
    if ( j > right_dist)
      low_row++;
    if (high_row < ulen)
      high_row ++;

    /* diagonalband */
    for (; i <= high_row; i++)
    {
      /* compute A_affine(i,j,I) */
      r_dist=add_safe_max(Atabcolumn[i][j-1].Rvalue,gap_extension+gap_opening);
      d_dist=add_safe_max(Atabcolumn[i][j-1].Dvalue,gap_extension+gap_opening);
      i_dist=add_safe_max(Atabcolumn[i][j-1].Ivalue,gap_extension);
      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[i][j].Ivalue = minvalue;
      Atabcolumn[i][j].Iedge = gt_linearalign_affinegapcost_set_edge(
                                                        r_dist, d_dist, i_dist);

      /* compute A_affine(i,j,R) */
      rcost = gt_scorehandler_get_replacement(scorehandler,
                                            useq[ustart+i-1], vseq[vstart+j-1]);
      r_dist = add_safe_max(Atabcolumn[i-1][j-1].Rvalue, rcost);
      d_dist = add_safe_max(Atabcolumn[i-1][j-1].Dvalue, rcost);
      i_dist = add_safe_max(Atabcolumn[i-1][j-1].Ivalue, rcost);
      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[i][j].Rvalue = minvalue;
      Atabcolumn[i][j].Redge = gt_linearalign_affinegapcost_set_edge(
                                                        r_dist, d_dist, i_dist);

      /* compute A_affine(i,j,D) */
      r_dist = add_safe_max(Atabcolumn[i-1][j].Rvalue,
                         gap_extension+gap_opening);
      d_dist = add_safe_max(Atabcolumn[i-1][j].Dvalue,gap_extension);
      i_dist = add_safe_max(Atabcolumn[i-1][j].Ivalue,
                          gap_extension+gap_opening);
      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[i][j].Dvalue = minvalue;
      Atabcolumn[i][j].Dedge = gt_linearalign_affinegapcost_set_edge(
                                                        r_dist, d_dist, i_dist);

    }
    /* above diagonal band */
    for (; i <= ulen; i++)
    {
      Atabcolumn[i][j].Rvalue = GT_WORD_MAX;
      Atabcolumn[i][j].Dvalue = GT_WORD_MAX;
      Atabcolumn[i][j].Ivalue = GT_WORD_MAX;
    }
  }
}

/* calculate alignment with diagonalband in square space with
 * affine gapcosts */
GtWord gt_diagonalbandalign_affinegapcost_in_square_space_generic(
                                             GtLinspaceManagement *space,
                                             const GtScoreHandler *scorehandler,
                                             GtAlignment *align,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen,
                                             GtWord left_dist,
                                             GtWord right_dist)
{
  GtWord distance;
  GtUword idx;
  GtAffinealignDPentry **Atabcolumn;

  gt_assert(align && scorehandler);
  if (space == NULL)
  {
    gt_array2dim_malloc(Atabcolumn, (ulen+1), (vlen+1));
  }
  else
  {
    gt_assert((ulen+1)*(vlen+1)*sizeof(**Atabcolumn) <=
               gt_linspace_management_get_valueTabsize(space));

    Atabcolumn = gt_linspace_management_get_rTabspace(space);
    *Atabcolumn = gt_linspace_management_get_valueTabspace(space);

    for (idx=1; idx<ulen+1; idx++)
      Atabcolumn[idx]=Atabcolumn[idx-1]+vlen+1;
  }

  diagonalband_fillDPtab_affine(Atabcolumn, useq, ustart, ulen, vseq, vstart,
                                vlen, left_dist, right_dist,
                                Affine_X, Affine_X, scorehandler);

  distance = MIN3(Atabcolumn[ulen][vlen].Rvalue,
                  Atabcolumn[ulen][vlen].Dvalue,
                  Atabcolumn[ulen][vlen].Ivalue);

  /* reconstruct alignment from 2dimarray Atabcolumn */
  gt_affinealign_traceback(align, Atabcolumn, ulen, vlen);

  if (space == NULL)
  {
    gt_array2dim_delete(Atabcolumn);
  }
  return distance;
}

/* calculate alignment with diagonalband in square space with
 * affine gapcosts */
GtWord gt_diagonalbandalign_affinegapcost_in_square_space(
                                                    GtLinspaceManagement *space,
                                                    GtAlignment *align,
                                                    const GtUchar *useq,
                                                    GtUword ustart,
                                                    GtUword ulen,
                                                    const GtUchar *vseq,
                                                    GtUword vstart,
                                                    GtUword vlen,
                                                    GtWord left_dist,
                                                    GtWord right_dist,
                                                    GtUword matchcost,
                                                    GtUword mismatchcost,
                                                    GtUword gap_opening,
                                                    GtUword gap_extension)
{
  GtWord distance;
  GtScoreHandler *scorehandler = gt_scorehandler_new(matchcost,mismatchcost,
                                                    gap_opening, gap_extension);

  distance = gt_diagonalbandalign_affinegapcost_in_square_space_generic(space,
                                     scorehandler, align, useq, ustart, ulen,
                                     vseq, vstart, vlen, left_dist, right_dist);
  gt_scorehandler_delete(scorehandler);

  return distance;
}

/* calculate only distance with diagonalband in square space  with
 * affine gapcosts */
GtWord gt_diagonalbandalign_affinegapcost_square_space_distance_only(
                                                           const GtUchar *useq,
                                                           GtUword ustart,
                                                           GtUword ulen,
                                                           const GtUchar *vseq,
                                                           GtUword vstart,
                                                           GtUword vlen,
                                                           GtWord left_dist,
                                                           GtWord right_dist,
                                                           const GtScoreHandler
                                                           *scorehandler)
{
  GtUword  distance;
  GtAffinealignDPentry **Atabcolumn;
  gt_assert(scorehandler);

   if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    return GT_WORD_MAX;
  }

  gt_array2dim_malloc(Atabcolumn, (ulen+1), (vlen+1));
  diagonalband_fillDPtab_affine(Atabcolumn, useq, ustart, ulen, vseq, vstart,
                                vlen, left_dist, right_dist,
                                Affine_X, Affine_X, scorehandler);

  distance = MIN3(Atabcolumn[ulen][vlen].Rvalue,
                  Atabcolumn[ulen][vlen].Dvalue,
                  Atabcolumn[ulen][vlen].Ivalue);

  gt_array2dim_delete(Atabcolumn);
  return distance;
}

static GtAffineAlignRnode evaluate_affineDBcrosspoints_from_2dimtab(
                                              GtAffineDiagAlignentry *Dtab,
                                              GtAffinealignDPentry **Atabcolumn,
                                              GtUword ulen, GtUword vlen,
                                              GtUword gap_opening,
                                              GtUword rowoffset,
                                              GtAffineAlignEdge from_edge,
                                              GtAffineAlignEdge edge)
{
  GtUword i, j;
  GtAffineAlignRnode rnode;
  GtDiagAlignentry *tempnode;
  gt_assert(Atabcolumn != NULL);

  i = ulen;
  j = vlen;

  edge = gt_linearalign_affinegapcost_minAdditionalCosts(&Atabcolumn[i][j],
                                                         edge, gap_opening);

  switch (edge)
  {
    case Affine_I:
      tempnode = &Dtab[vlen].val_I;
      rnode = (GtAffineAlignRnode) {vlen, Affine_I};
      break;
    case Affine_D:
      tempnode = &Dtab[vlen].val_D;
      rnode = (GtAffineAlignRnode) {vlen, Affine_D};
      break;
    default:
      tempnode = &Dtab[vlen].val_R;
      rnode = (GtAffineAlignRnode) {vlen, Affine_R};
  }

  while (i > 0 || j > 0) {
    if (j == vlen)
      rnode.edge = edge;
    switch (edge) {
      case Affine_R:
        gt_assert(Atabcolumn[i][j].Rvalue != GT_WORD_MAX);
        Dtab[j].val_R.currentrowindex = i + rowoffset;
        edge = Atabcolumn[i][j].Redge;
        tempnode->last_type = Affine_R;
        tempnode = &Dtab[j].val_R;
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
        Dtab[j].val_I.currentrowindex = i + rowoffset;
        edge = Atabcolumn[i][j].Iedge;
        tempnode->last_type = Affine_I;
        tempnode = &Dtab[j].val_I;
        gt_assert(j);
        j--;
        break;
      default:
        gt_assert(false);
    }
  }
  tempnode->last_type = edge;
  /* special case for first crosspoint */
  Dtab[0].val_R = (GtDiagAlignentry) {GT_UWORD_MAX, rowoffset, from_edge};
  Dtab[0].val_D = (GtDiagAlignentry) {GT_UWORD_MAX, rowoffset, from_edge};
  Dtab[0].val_I = (GtDiagAlignentry) {GT_UWORD_MAX, rowoffset, from_edge};

  return rnode;
}

/* create affine DBcrosspointtab to combine square calculating with linear
 * calculating. from_edge describes type of crosspoint node, edge describes the
 * incoming way to next unkonown crosspoint and to_edge describes type of
 * previous crosspoint.
 * Returns edge and index of lastcrosspoint in matrix.
 */
static GtAffineAlignRnode affineDtab_in_square_space(
                                                   GtLinspaceManagement *space,
                                                   GtAffineDiagAlignentry *Dtab,
                                                   const GtUchar *useq,
                                                   GtUword ustart,
                                                   GtUword ulen,
                                                   const GtUchar *vseq,
                                                   GtUword vstart,
                                                   GtUword vlen,
                                                   GtWord left_dist,
                                                   GtWord right_dist,
                                                   GtUword rowoffset,
                                                   GtAffineAlignEdge from_edge,
                                                   GtAffineAlignEdge edge,
                                                   GtAffineAlignEdge to_edge,
                                                   const GtScoreHandler
                                                   *scorehandler)
{
  GtAffinealignDPentry **Atabcolumn;
  GtUword idx, gap_opening;

  gt_assert(Dtab && space && scorehandler);

  gt_assert((ulen+1)*(vlen+1)*sizeof(**Atabcolumn) <=
               gt_linspace_management_get_valueTabsize(space));

  Atabcolumn = gt_linspace_management_get_rTabspace(space);
  *Atabcolumn = gt_linspace_management_get_valueTabspace(space);

  for (idx=1;idx<ulen+1;idx++)
    Atabcolumn[idx]=Atabcolumn[idx-1]+vlen+1;
  diagonalband_fillDPtab_affine(Atabcolumn, useq, ustart, ulen, vseq, vstart,
                                vlen, left_dist, right_dist,
                                from_edge, edge, scorehandler);

  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  GtAffineAlignRnode rnode = evaluate_affineDBcrosspoints_from_2dimtab(Dtab,
                                            Atabcolumn, ulen, vlen, gap_opening,
                                            rowoffset, from_edge, to_edge);

  return rnode;
}

/* calculate only distance with diagonalband in linear space O(n)
 * with affine gapcosts */
static GtWord diagonalband_linear_affine(const GtUchar *useq,
                                          GtUword ustart,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vstart,
                                          GtUword vlen,
                                          GtWord left_dist,
                                          GtWord right_dist,
                                          const GtScoreHandler *scorehandler)
{
  GtUword colindex, rowindex, low_row, high_row, width,
          gap_opening, gap_extension;
  GtWord distance, rcost, r_dist, d_dist, i_dist, minvalue;
  GtAffinealignDPentry *Atabcolumn, northwestAffinealignDPentry,
                       westAffinealignDPentry;
  bool last_row = false;
  distance = GT_WORD_MAX;

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }
  gt_assert(scorehandler);
  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);

  width = right_dist - left_dist + 1;
  Atabcolumn = gt_malloc(sizeof(*Atabcolumn) * width);

  low_row = 0;
  high_row = -left_dist;
  Atabcolumn[low_row].Rvalue = 0;
  Atabcolumn[low_row].Dvalue = gap_opening;
  Atabcolumn[low_row].Ivalue = gap_opening;

  for (rowindex = low_row+1; rowindex <= high_row; rowindex ++)
  {
    Atabcolumn[rowindex-low_row].Rvalue = GT_WORD_MAX;
    Atabcolumn[rowindex-low_row].Dvalue = add_safe_max(
                                          Atabcolumn[rowindex-low_row-1].Dvalue,
                                          gap_extension);
    Atabcolumn[rowindex-low_row].Ivalue = GT_WORD_MAX;
  }
  if (high_row == ulen)
    last_row = true;
  for (colindex = 1; colindex <= vlen; colindex++)
  {
    northwestAffinealignDPentry = Atabcolumn[0];
    if (colindex > right_dist)
    {
      if (low_row != high_row)
        westAffinealignDPentry = Atabcolumn[1];
      low_row++;
    }
    else
      westAffinealignDPentry = Atabcolumn[0];

    if (high_row < ulen)
      high_row ++;
    if (!last_row && rowindex == high_row)
      {
        westAffinealignDPentry.Rvalue = GT_WORD_MAX;
        westAffinealignDPentry.Dvalue = GT_WORD_MAX;
        westAffinealignDPentry.Ivalue = GT_WORD_MAX;
      }

    r_dist = add_safe_max(westAffinealignDPentry.Rvalue,
                          gap_extension+gap_opening);
    d_dist = add_safe_max(westAffinealignDPentry.Dvalue,
                          gap_extension+gap_opening);
    i_dist = add_safe_max(westAffinealignDPentry.Ivalue,gap_extension);
    minvalue = MIN3(r_dist, d_dist, i_dist);
    Atabcolumn[0].Ivalue = minvalue;
    Atabcolumn[0].Rvalue = GT_WORD_MAX;
    Atabcolumn[0].Dvalue = GT_WORD_MAX;

    if (low_row > 0 )
    {
      rcost = gt_scorehandler_get_replacement(scorehandler,
                              useq[ustart+rowindex-1],vseq[vstart+colindex-1]);

      r_dist = add_safe_max(northwestAffinealignDPentry.Rvalue, rcost);
      d_dist = add_safe_max(northwestAffinealignDPentry.Dvalue, rcost);
      i_dist = add_safe_max(northwestAffinealignDPentry.Ivalue, rcost);
      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[0].Rvalue = minvalue;
    }
    for (rowindex = low_row + 1; rowindex <= high_row; rowindex++)
    {
      northwestAffinealignDPentry = westAffinealignDPentry;
      if (!last_row && rowindex == high_row)
      {
        westAffinealignDPentry.Rvalue = GT_WORD_MAX;
        westAffinealignDPentry.Dvalue = GT_WORD_MAX;
        westAffinealignDPentry.Ivalue = GT_WORD_MAX;
      }
      else if (low_row > 0)
        westAffinealignDPentry = Atabcolumn[rowindex-low_row+1];
      else
        westAffinealignDPentry = Atabcolumn[rowindex-low_row];

      if (rowindex == ulen)
        last_row = true;

      r_dist = add_safe_max(westAffinealignDPentry.Rvalue,
                            gap_extension+gap_opening);
      d_dist = add_safe_max(westAffinealignDPentry.Dvalue,
                            gap_extension+gap_opening);
      i_dist = add_safe_max(westAffinealignDPentry.Ivalue,gap_extension);

      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[rowindex-low_row].Ivalue = minvalue;

      rcost = gt_scorehandler_get_replacement(scorehandler,
                              useq[ustart+rowindex-1], vseq[vstart+colindex-1]);

      r_dist = add_safe_max(northwestAffinealignDPentry.Rvalue, rcost);
      d_dist = add_safe_max(northwestAffinealignDPentry.Dvalue, rcost);
      i_dist = add_safe_max(northwestAffinealignDPentry.Ivalue, rcost);

      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[rowindex-low_row].Rvalue = minvalue;

      r_dist = add_safe_max(Atabcolumn[rowindex-low_row-1].Rvalue,
                            gap_extension+gap_opening);
      d_dist = add_safe_max(Atabcolumn[rowindex-low_row-1].Dvalue,
                           gap_extension);
      i_dist = add_safe_max(Atabcolumn[rowindex-low_row-1].Ivalue,
                           gap_extension+gap_opening);

      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[rowindex-low_row].Dvalue = minvalue;
    }
  }

  distance = MIN3(Atabcolumn[high_row-low_row].Rvalue,
                  Atabcolumn[high_row-low_row].Dvalue,
                  Atabcolumn[high_row-low_row].Ivalue);
  gt_free(Atabcolumn);

  return distance;
}

/* helpfunctions */
static void inline  set_invalid_Diagentry(GtDiagAlignentry *node)
{
  gt_assert(node != NULL);
  node->currentrowindex = GT_UWORD_MAX;
  node->last_type = Affine_X;
  node->lastcpoint = GT_UWORD_MAX;
}

static void inline set_valid_Diagentry(GtDiagAlignentry *node_to,
                                       const GtAffineAlignRtabentry *entry_from,
                                       GtWord minvalue, GtWord r_dist,
                                       GtWord i_dist, GtWord d_dist)
{
  gt_assert(node_to != NULL && entry_from != NULL);
  if (minvalue == r_dist)
  {
    node_to->last_type = entry_from->val_R.edge;
    node_to->lastcpoint = entry_from->val_R.idx;
  }
  else if (minvalue == i_dist)
  {
    node_to->last_type = entry_from->val_I.edge;
    node_to->lastcpoint = entry_from->val_I.idx;
  }
  else if (minvalue == d_dist)
  {
    node_to->last_type = entry_from->val_D.edge;
    node_to->lastcpoint = entry_from->val_D.idx;
  }

}

static void inline set_invalid_Rnode(GtAffineAlignRnode *node)
{
  gt_assert(node != NULL);
  node->idx = GT_UWORD_MAX;
  node->edge = Affine_X;
}

static void inline set_valid_Rnode(GtAffineAlignRnode *node_to,
                                   GtAffineAlignRtabentry *entry_from,
                                   GtWord minvalue, GtWord r_dist,
                                   GtWord i_dist, GtWord d_dist)
{
  gt_assert(node_to != NULL && entry_from != NULL);
  if (minvalue == r_dist)
    *node_to = entry_from->val_R;
  else if (minvalue == i_dist)
     *node_to = entry_from->val_I;
  else if (minvalue == d_dist)
     *node_to = entry_from->val_D;
}

/* calculate first column */
static void firstaffineDBtabcolumn(GtAffinealignDPentry *Atabcolumn,
                                   GtAffineAlignRtabentry *Rtabcolumn,
                                   GtAffineDiagAlignentry *Diagcolumn,
                                   GtAffineAlignEdge edge,
                                   GtAffineAlignEdge from_edge,
                                   GtUword offset,
                                   GtWord left_dist,
                                   GtWord right_dist,
                                   GtUword gap_opening,
                                   GtUword gap_extension)
{
  GtUword rowindex, low_row, high_row;
  GtWord diag;

  diag = GT_DIV2(left_dist + right_dist);
  low_row = 0;
  high_row = -left_dist;

  Atabcolumn[low_row].Rvalue = GT_WORD_MAX;
  Atabcolumn[low_row].Dvalue = GT_WORD_MAX;
  Atabcolumn[low_row].Ivalue = GT_WORD_MAX;

  set_invalid_Diagentry(&Diagcolumn[0].val_R);
  set_invalid_Diagentry(&Diagcolumn[0].val_D);
  set_invalid_Diagentry(&Diagcolumn[0].val_I);

  set_invalid_Rnode(&Rtabcolumn[0].val_R);
  set_invalid_Rnode(&Rtabcolumn[0].val_D);
  set_invalid_Rnode(&Rtabcolumn[0].val_I);

  switch (edge) {
  case Affine_R:
    Atabcolumn[low_row].Rvalue = 0;
    Rtabcolumn[0].val_R.edge = from_edge;
    if (diag == 0)
    {
      Diagcolumn[0].val_R.currentrowindex = 0 + offset;
      Diagcolumn[0].val_R.last_type = from_edge;
      Rtabcolumn[0].val_R.idx = 0;
      Rtabcolumn[0].val_R.edge = Affine_R;
    }
    break;
  case Affine_D:
    Atabcolumn[low_row].Dvalue = 0;
    Rtabcolumn[0].val_D.edge = from_edge;
    if (diag == 0)
    {
      Diagcolumn[0].val_D.currentrowindex = 0 + offset;
      Diagcolumn[0].val_D.last_type = from_edge;
      Rtabcolumn[0].val_D.idx = 0;
      Rtabcolumn[0].val_D.edge = Affine_D;
    }
    break;
  case Affine_I:
    Atabcolumn[low_row].Ivalue = 0;
    Rtabcolumn[0].val_I.edge = from_edge;
    if (diag == 0)
    {
      Diagcolumn[0].val_I.currentrowindex = 0 + offset;
      Diagcolumn[0].val_I.last_type = from_edge;
      Rtabcolumn[0].val_I.idx = 0;
      Rtabcolumn[0].val_I.edge = Affine_I;
    }
    break;
  default:
    Atabcolumn[low_row].Rvalue = 0;
    Atabcolumn[low_row].Dvalue = gap_opening;
    Atabcolumn[low_row].Ivalue = gap_opening;
    Rtabcolumn[0].val_I.edge = from_edge;
    Rtabcolumn[0].val_R.edge = from_edge;
    Rtabcolumn[0].val_D.edge = from_edge;
    if (diag == 0)
    {
      Diagcolumn[0].val_R.currentrowindex = 0 + offset;
      Diagcolumn[0].val_D.currentrowindex = 0 + offset;
      Diagcolumn[0].val_I.currentrowindex = 0 + offset;

      Rtabcolumn[0].val_R.idx = 0;
      Rtabcolumn[0].val_R.edge = Affine_R;
      Rtabcolumn[0].val_D.idx = 0;
      Rtabcolumn[0].val_D.edge = Affine_D;
      Rtabcolumn[0].val_I.idx = 0;
      Rtabcolumn[0].val_I.edge = Affine_I;
    }
  }

  for (rowindex = low_row+1; rowindex <= high_row; rowindex++)
  {
    Atabcolumn[rowindex-low_row].Rvalue = GT_WORD_MAX;
    Atabcolumn[rowindex-low_row].Dvalue = add_safe_max(
                                          Atabcolumn[rowindex-low_row-1].Dvalue,
                                          gap_extension);
    Atabcolumn[rowindex-low_row].Ivalue = GT_WORD_MAX;

    if (diag == -(GtWord)rowindex)
    {
      Diagcolumn[0].val_D.last_type = from_edge;
      Diagcolumn[0].val_D.lastcpoint = GT_UWORD_MAX;
      Diagcolumn[0].val_D.currentrowindex = rowindex + offset;
      Rtabcolumn[rowindex-low_row].val_D.idx = 0;
      Rtabcolumn[rowindex-low_row].val_D.edge = Affine_D;
      set_invalid_Diagentry(&Diagcolumn[0].val_R);
      set_invalid_Diagentry(&Diagcolumn[0].val_I);
    }
    else
    {
      Rtabcolumn[rowindex-low_row] = Rtabcolumn[rowindex-low_row-1];
    }
  }

}

/* calculate all columns */
static GtAffineAlignRnode evaluateallaffineDBcolumns(
                                            GtLinspaceManagement *spacemanager,
                                            GtAffineDiagAlignentry *Diagcolumn,
                                            const GtScoreHandler *scorehandler,
                                            GtAffineAlignEdge edge,
                                            GtAffineAlignEdge from_edge,
                                            GtAffineAlignEdge to_edge,
                                            GtUword offset,
                                            const GtUchar *useq,
                                            GtUword ustart, GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart, GtUword vlen,
                                            GtWord left_dist, GtWord right_dist)
{
  GtUword gap_extension, gap_opening, colindex, rowindex, low_row, high_row;
  /*lowest and highest row between a diagonal band*/
  GtWord diag, r_dist, d_dist, i_dist, minvalue, rcost;

  bool last_row = false;
  GtAffinealignDPentry *Atabcolumn, northwestAffinealignDPentry,
  westAffinealignDPentry = (GtAffinealignDPentry)
                           {GT_WORD_MAX, GT_WORD_MAX, GT_WORD_MAX};;
  GtAffineAlignRtabentry *Rtabcolumn, northwestRtabentry, westRtabentry = {{0}};
  GtAffineAlignRnode lastcpoint = {GT_UWORD_MAX, Affine_X};

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }

  Atabcolumn = gt_linspace_management_get_valueTabspace(spacemanager);
  Rtabcolumn = gt_linspace_management_get_rTabspace(spacemanager);

  diag = GT_DIV2(left_dist + right_dist);
  low_row = 0;
  high_row = -left_dist;
  if (high_row == ulen)
    last_row = true;

  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);

 /* first column */
  firstaffineDBtabcolumn(Atabcolumn, Rtabcolumn, Diagcolumn, edge, from_edge,
                   offset, left_dist, right_dist, gap_opening, gap_extension);
  /* next columns */
  for (colindex = 1; colindex <= vlen; colindex++)
  {
    northwestAffinealignDPentry = Atabcolumn[0];
    northwestRtabentry = Rtabcolumn[0];

    if (colindex > right_dist)
    {
      if (low_row != high_row) {
        westAffinealignDPentry = Atabcolumn[1];
        westRtabentry = Rtabcolumn[1];
      }
      low_row++;
    }
    else
    {
      westAffinealignDPentry = Atabcolumn[0];
      westRtabentry = Rtabcolumn[0];
    }
    if (high_row < ulen)
      high_row ++;
    if (!last_row && low_row == high_row)
    {/* prev is outside of diagonalband*/
      westAffinealignDPentry = (GtAffinealignDPentry)
                               {GT_WORD_MAX, GT_WORD_MAX, GT_WORD_MAX};
      westRtabentry.val_R = (GtAffineAlignRnode) {GT_UWORD_MAX, Affine_X};
      westRtabentry.val_D = (GtAffineAlignRnode) {GT_UWORD_MAX, Affine_X};
      westRtabentry.val_I = (GtAffineAlignRnode) {GT_UWORD_MAX, Affine_X};
    }

    /* insertion */
    r_dist = add_safe_max(westAffinealignDPentry.Rvalue,
                          gap_extension + gap_opening);
    d_dist = add_safe_max(westAffinealignDPentry.Dvalue,
                          gap_extension + gap_opening);
    i_dist = add_safe_max(westAffinealignDPentry.Ivalue, gap_extension);

    minvalue = MIN3(r_dist, d_dist, i_dist);
    Atabcolumn[0].Ivalue = minvalue;
    Atabcolumn[0].Rvalue = GT_WORD_MAX;
    Atabcolumn[0].Dvalue = GT_WORD_MAX;

    if (diag == (GtWord)colindex - (GtWord)low_row)
    {
      set_invalid_Diagentry(&Diagcolumn[colindex].val_R);
      set_invalid_Diagentry(&Diagcolumn[colindex].val_D);
      set_valid_Diagentry(&Diagcolumn[colindex].val_I, &westRtabentry, minvalue,
                          r_dist, i_dist, d_dist);
      Diagcolumn[colindex].val_I.currentrowindex = low_row + offset;
      set_invalid_Rnode(&Rtabcolumn[0].val_R);
      set_invalid_Rnode(&Rtabcolumn[0].val_D);
      Rtabcolumn[0].val_I.idx = colindex;
      Rtabcolumn[0].val_I.edge = Affine_I;
    }
    else
    {
      set_valid_Rnode(&Rtabcolumn[0].val_I, &westRtabentry, minvalue,
                       r_dist, i_dist, d_dist);
      Rtabcolumn[0].val_D = (GtAffineAlignRnode) {GT_UWORD_MAX, Affine_X};
      Rtabcolumn[0].val_R = (GtAffineAlignRnode) {GT_UWORD_MAX, Affine_X};
    }

    /* replacement possible for 0-entry */
    if (low_row > 0 )
    {
      rcost = gt_scorehandler_get_replacement(scorehandler,
                               useq[ustart+low_row-1], vseq[vstart+colindex-1]);

      r_dist = add_safe_max(northwestAffinealignDPentry.Rvalue, rcost);
      d_dist = add_safe_max(northwestAffinealignDPentry.Dvalue, rcost);
      i_dist = add_safe_max(northwestAffinealignDPentry.Ivalue, rcost);

      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[0].Rvalue = minvalue;

      if (diag == (GtWord)colindex - (GtWord)low_row)
      {
        set_valid_Diagentry(&Diagcolumn[colindex].val_R, &northwestRtabentry,
                            minvalue, r_dist, i_dist, d_dist);
        Diagcolumn[colindex].val_R.currentrowindex = low_row + offset;
        Rtabcolumn[0].val_R.idx = colindex;
        Rtabcolumn[0].val_R.edge = Affine_R;
      }
      else
      {
        set_valid_Rnode(&Rtabcolumn[0].val_R, &northwestRtabentry,
                         minvalue, r_dist, i_dist, d_dist);
      }
    }
    for (rowindex = low_row + 1; rowindex <= high_row; rowindex++)
    {
      northwestAffinealignDPentry = westAffinealignDPentry;
      northwestRtabentry = westRtabentry;

      if (!last_row && rowindex == high_row)
      {/* prev is outside of diagonalband*/
        westAffinealignDPentry = (GtAffinealignDPentry)
                                 {GT_WORD_MAX, GT_WORD_MAX, GT_WORD_MAX};
        westRtabentry.val_R = (GtAffineAlignRnode) {GT_UWORD_MAX, Affine_X};
        westRtabentry.val_D = (GtAffineAlignRnode) {GT_UWORD_MAX, Affine_X};
        westRtabentry.val_I = (GtAffineAlignRnode) {GT_UWORD_MAX, Affine_X};
      }
      else if (low_row > 0)
      {/* shifted diagonalband*/
        westAffinealignDPentry = Atabcolumn[rowindex-low_row+1];
        westRtabentry = Rtabcolumn[rowindex-low_row+1];
      }
      else
      {/* normaly prev*/
        westAffinealignDPentry = Atabcolumn[rowindex-low_row];
        westRtabentry = Rtabcolumn[rowindex-low_row];
      }
      if (rowindex == ulen)
        last_row = true;

      /* insertion */
      r_dist = add_safe_max(westAffinealignDPentry.Rvalue,
                            gap_extension+gap_opening);
      d_dist = add_safe_max(westAffinealignDPentry.Dvalue,
                            gap_extension+gap_opening);
      i_dist = add_safe_max(westAffinealignDPentry.Ivalue,gap_extension);

      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[rowindex-low_row].Ivalue = minvalue;
      if (diag == (GtWord)colindex - (GtWord)rowindex)
      {
        set_valid_Diagentry(&Diagcolumn[colindex].val_I, &westRtabentry,
                             minvalue, r_dist, i_dist, d_dist);
        Diagcolumn[colindex].val_I.currentrowindex = rowindex+offset;

        Rtabcolumn[rowindex-low_row].val_I.idx = colindex;
        Rtabcolumn[rowindex-low_row].val_I.edge = Affine_I;
      }
      else
      {
        set_valid_Rnode(&Rtabcolumn[rowindex-low_row].val_I, &westRtabentry,
                         minvalue,r_dist,i_dist,d_dist);
      }
      /* replacement */
      rcost = gt_scorehandler_get_replacement(scorehandler,
                             useq[ustart+rowindex-1], vseq[vstart+colindex-1]);

      r_dist = add_safe_max(northwestAffinealignDPentry.Rvalue, rcost);
      d_dist = add_safe_max(northwestAffinealignDPentry.Dvalue, rcost);
      i_dist = add_safe_max(northwestAffinealignDPentry.Ivalue, rcost);
      minvalue = MIN3(r_dist, d_dist, i_dist);

      Atabcolumn[rowindex-low_row].Rvalue = minvalue;
      if (diag == (GtWord)colindex - (GtWord)rowindex)
      {
        set_valid_Diagentry(&Diagcolumn[colindex].val_R, &northwestRtabentry,
                            minvalue, r_dist, i_dist, d_dist);
        Diagcolumn[colindex].val_R.currentrowindex = rowindex+offset;
        Rtabcolumn[rowindex-low_row].val_R.idx = colindex;
        Rtabcolumn[rowindex-low_row].val_R.edge = Affine_R;
      }
      else
      {
        set_valid_Rnode(&Rtabcolumn[rowindex-low_row].val_R,&northwestRtabentry,
                         minvalue,r_dist,i_dist,d_dist);
      }

      /* deletion */
      r_dist = add_safe_max(Atabcolumn[rowindex-low_row-1].Rvalue,
                           gap_extension+gap_opening);
      d_dist = add_safe_max(Atabcolumn[rowindex-low_row-1].Dvalue,
                           gap_extension);
      i_dist = add_safe_max(Atabcolumn[rowindex-low_row-1].Ivalue,
                           gap_extension+gap_opening);

      minvalue = MIN3(r_dist, d_dist, i_dist);
      Atabcolumn[rowindex-low_row].Dvalue = minvalue;

      if (diag == (GtWord)colindex - (GtWord)rowindex)
      {
        set_valid_Diagentry(&Diagcolumn[colindex].val_D,
                           &Rtabcolumn[rowindex-low_row-1], minvalue,
                            r_dist,i_dist, d_dist);
        Diagcolumn[colindex].val_D.currentrowindex = rowindex+offset;
        Rtabcolumn[rowindex-low_row].val_D.idx = colindex;
        Rtabcolumn[rowindex-low_row].val_D.edge = Affine_D;
      }
      else
      {
        set_valid_Rnode(&Rtabcolumn[rowindex-low_row].val_D,
                        &Rtabcolumn[rowindex-low_row-1], minvalue,
                         r_dist,i_dist,d_dist);
      }
    }

  }
  /* last crosspoint of optimal path */
  r_dist = Atabcolumn[high_row-low_row].Rvalue;
  d_dist = Atabcolumn[high_row-low_row].Dvalue;
  i_dist = Atabcolumn[high_row-low_row].Ivalue;

  switch (to_edge)
  {
    case Affine_I:
      r_dist = add_safe_max (r_dist, gap_opening);
      d_dist = add_safe_max (d_dist, gap_opening);
      break;
    case Affine_D:
      r_dist = add_safe_max (r_dist, gap_opening);
      i_dist = add_safe_max (i_dist, gap_opening);
      break;
    default:
      break;
  }

  minvalue = MIN3(r_dist, d_dist, i_dist);
  if (minvalue == r_dist)
    lastcpoint = Rtabcolumn[high_row-low_row].val_R;
  else if (minvalue == i_dist)
    lastcpoint = Rtabcolumn[high_row-low_row].val_I;
  else if (minvalue == d_dist)
    lastcpoint = Rtabcolumn[high_row-low_row].val_D;

  return lastcpoint;
}

/* calculate affine crosspoint realting to diagonal in recursive way */
static GtAffineAlignRnode evaluateaffineDBcrosspoints(
                                             GtLinspaceManagement *spacemanager,
                                             GtAffineDiagAlignentry *Diagcolumn,
                                             const GtScoreHandler *scorehandler,
                                             GtAffineAlignEdge edge,
                                             GtAffineAlignEdge from_edge,
                                             GtAffineAlignEdge to_edge,
                                             GtUword rowoffset,
                                             GtUword coloffset,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen,
                                             GtWord left_dist,
                                             GtWord right_dist)
{
  GtUword i,new_ulen, new_vlen, col_start, col_end,
          row_start = 0, row_end;
  GtWord new_left, new_right, diag;
  GtDiagAlignentry cpoint = {0,0,0}, prevcpoint;
  GtAffineDiagAlignentry temp_entry;
  GtAffineAlignRnode rpoint, temprpoint, lastrpoint;
  GtAffineAlignEdge prevcp_type,cp_type;

  diag = GT_DIV2(left_dist+right_dist);
  gt_assert(vstart == coloffset);

  if (ulen == 0)
  {
    switch (edge) {
      case(Affine_R):
      {
        Diagcolumn[0].val_R.currentrowindex = rowoffset;
        Diagcolumn[0].val_R.last_type = from_edge;
        prevcp_type = Affine_R;
        break;
      }
      case(Affine_I):
      {
        Diagcolumn[0].val_I.currentrowindex = rowoffset;
        Diagcolumn[0].val_I.last_type = from_edge;
        prevcp_type = Affine_I;
        break;
      }
      case(Affine_D):
      {
        Diagcolumn[0].val_D.currentrowindex = rowoffset;
        Diagcolumn[0].val_D.last_type = from_edge;
        prevcp_type = Affine_D;
        break;
      }
      default:
      {
        Diagcolumn[0].val_I.currentrowindex = rowoffset;
        Diagcolumn[0].val_I.last_type = from_edge;
        prevcp_type = Affine_I;
      }
    }

    for (i = 1; i <=vlen; i++)
    {
      Diagcolumn[i].val_I.currentrowindex = rowoffset;
      Diagcolumn[i].val_I.last_type = prevcp_type;
      prevcp_type = Affine_I;
    }
    return (GtAffineAlignRnode) {vlen, prevcp_type};

  }

  if (vlen == 0)
  {
    switch (edge) {
      case(Affine_R):
      {
        Diagcolumn[0].val_R.currentrowindex = rowoffset;
        Diagcolumn[0].val_R.last_type = from_edge;
        break;
      }
      case(Affine_I):
      {
        Diagcolumn[0].val_I.currentrowindex = rowoffset;
        Diagcolumn[0].val_I.last_type = from_edge;
        break;
      }
      case(Affine_D):
      {
        Diagcolumn[0].val_D.currentrowindex = rowoffset;
        Diagcolumn[0].val_D.last_type = from_edge;
        break;
      }
      default:
      {
        Diagcolumn[0].val_D.currentrowindex = rowoffset;
        Diagcolumn[0].val_D.last_type = from_edge;
      }
    }
    return (GtAffineAlignRnode) {0, edge};
  }

  if (gt_linspace_management_checksquare(spacemanager, ulen, vlen,
                                         sizeof (GtAffinealignDPentry),
                                         sizeof (GtAffineAlignRtabentry)))
  { /* call square function */
    return affineDtab_in_square_space(spacemanager, Diagcolumn,
                                      useq, ustart, ulen, vseq, vstart, vlen,
                                      left_dist, right_dist, rowoffset,
                                      from_edge, edge, to_edge,scorehandler);
  }

  rpoint = evaluateallaffineDBcolumns(spacemanager, Diagcolumn, scorehandler,
                                      edge, from_edge, to_edge, rowoffset,
                                      useq, ustart, ulen,  vseq, vstart, vlen,
                                      left_dist, right_dist);

  lastrpoint = rpoint;
  col_start = rpoint.idx;
  cp_type = rpoint.edge;

  /* if no crosspoint is found */
  if (col_start == GT_UWORD_MAX)
  {
    gt_assert(diag != 0);
    if (diag < 0)
    {
      return evaluateaffineDBcrosspoints(spacemanager, Diagcolumn, scorehandler,
                                         edge, from_edge, to_edge, rowoffset,
                                         coloffset, useq, ustart, ulen, vseq,
                                         vstart, vlen, diag+1, right_dist);
    } else
    {
      if (diag > 0)
      {
        return evaluateaffineDBcrosspoints(spacemanager, Diagcolumn,
                                           scorehandler,
                                           edge, from_edge, to_edge, rowoffset,
                                           coloffset, useq, ustart, ulen, vseq,
                                           vstart, vlen, left_dist, diag-1);
      }
    }
  }
  else
  {
    switch (cp_type) {
    case Affine_R:
      cpoint = Diagcolumn[col_start].val_R;
      cp_type = Affine_R;
      row_start = Diagcolumn[col_start].val_R.currentrowindex;
      break;
    case Affine_D:
      cpoint = Diagcolumn[col_start].val_D;
      cp_type = Affine_D;
      row_start = Diagcolumn[col_start].val_D.currentrowindex;
      break;
    default:
      gt_assert(cp_type == Affine_I);
      cpoint = Diagcolumn[col_start].val_I;
      cp_type = Affine_I;
      row_start = Diagcolumn[col_start].val_I.currentrowindex;
      break;
    }
  }

  /* exception, if last cpoint != (m+1)entry */
  if (col_start != vlen)
  {
    if (diag + ((GtWord)ulen-(GtWord)vlen) > 0)
    {
      new_ulen = ulen - (row_start+1-rowoffset);
      new_vlen = vlen-col_start;
      new_left = MAX((GtWord)left_dist-diag+1,
                    -(GtWord)new_ulen);
      new_right = 0;
      temp_entry = Diagcolumn[col_start];
      lastrpoint = evaluateaffineDBcrosspoints(spacemanager,
                      Diagcolumn+col_start, scorehandler, Affine_D,
                      cpoint.last_type, to_edge,
                      row_start+1, coloffset+col_start, useq,
                      row_start+1, new_ulen, vseq,
                      vstart+col_start, new_vlen, new_left, new_right);
      Diagcolumn[col_start] = temp_entry;
      Diagcolumn[col_start+1].val_R.last_type = cp_type;
      Diagcolumn[col_start+1].val_D.last_type = cp_type;
      Diagcolumn[col_start+1].val_I.last_type = cp_type;

      lastrpoint.idx += col_start;
    }
    else
    {
      new_ulen = ulen - (row_start - rowoffset);
      new_vlen = vlen-col_start-1;
      new_left = -1;
      new_right = MIN((GtWord)right_dist-((GtWord)diag)-1,new_vlen);

      lastrpoint = evaluateaffineDBcrosspoints(spacemanager,
                            Diagcolumn+col_start+1,scorehandler, Affine_I,
                            cp_type, to_edge, row_start, coloffset+col_start+1,
                            useq, row_start, new_ulen, vseq, vstart+col_start+1,
                            new_vlen, new_left, new_right);
      lastrpoint.idx += (col_start+1);
    }
  }
  /* look at all 'normally' crosspoints */
  while (cpoint.lastcpoint != GT_UWORD_MAX)
  {
    prevcpoint = cpoint;
    prevcp_type = cp_type;
    col_end = col_start;
    row_end = row_start;
    col_start = prevcpoint.lastcpoint;
    switch (prevcpoint.last_type) {
    case Affine_R:
      cpoint = Diagcolumn[col_start].val_R;
      cp_type = Affine_R;
      row_start = Diagcolumn[col_start].val_R.currentrowindex;
      break;
    case Affine_D:
      cpoint = Diagcolumn[col_start].val_D;
      cp_type = Affine_D;
      row_start = Diagcolumn[col_start].val_D.currentrowindex;
      break;
    default:
      gt_assert(prevcpoint.last_type == Affine_I);
      cpoint = Diagcolumn[col_start].val_I;
      cp_type = Affine_I;
      row_start = Diagcolumn[col_start].val_I.currentrowindex;
      break;
    }

    if (prevcp_type == Affine_R ||
       ((prevcp_type == Affine_I) && (col_end-col_start == 1)))
    {
      continue;/* next crosspoint is also on the diagonal*/
    }
    else if (prevcp_type == Affine_D)
    {
      new_ulen = row_end - row_start - 1;
      new_vlen = col_end-col_start-1;
      new_left = -1;
      new_right = MIN(right_dist-diag-1,new_vlen);

      temprpoint = evaluateaffineDBcrosspoints(spacemanager,
                           Diagcolumn+col_start+1, scorehandler,
                           Affine_I, cp_type, Affine_D, row_start,
                           coloffset + col_start + 1, useq, row_start,
                           new_ulen, vseq, vstart + col_start + 1,
                           new_vlen, new_left, new_right);
      if (temprpoint.idx + col_start + 1 < vlen)
      {
        GtUword update_idx = temprpoint.idx+1+col_start+1;
        Diagcolumn[update_idx].val_R.last_type = temprpoint.edge;
        Diagcolumn[update_idx].val_D.last_type = temprpoint.edge;
        Diagcolumn[update_idx].val_I.last_type = temprpoint.edge;
      }
      if (temprpoint.idx + col_start + 1 == lastrpoint.idx)
      {
        lastrpoint = temprpoint;
        lastrpoint.idx += col_start + 1;
      }
    }
    else if (prevcp_type == Affine_I)
    {
      new_ulen = row_end - row_start-1;
      new_left = MAX(left_dist-diag+1,
                        -(GtWord)new_ulen);
      new_right = 0;
      temp_entry = Diagcolumn[col_start];
      temprpoint = evaluateaffineDBcrosspoints(spacemanager,
                            Diagcolumn + col_start, scorehandler,
                            Affine_D, cpoint.last_type, Affine_I, row_start + 1,
                            coloffset + col_start,
                            useq, row_start+1, new_ulen,
                            vseq, vstart+col_start, col_end-col_start-1,
                            new_left, new_right);
      Diagcolumn[col_start] = temp_entry;
      Diagcolumn[col_start+1].val_R.last_type = cp_type;
      Diagcolumn[col_start+1].val_D.last_type = cp_type;
      Diagcolumn[col_start+1].val_I.last_type = cp_type;
      Diagcolumn[col_end].val_I.last_type = temprpoint.edge;
    }
    else
    {
      /* if (Diagcolumn[cpoint].last_type == Linear_X), never reach this line */
      gt_assert(false);
    }
  }
  col_end = col_start;
  row_end = row_start;
  /* exception, if first crosspoint != 0-entry */
   if (vstart-coloffset != col_end)
   {
     switch (cp_type) {
     case Affine_D:
      if (row_end == ustart-1)
         gt_assert(false);
       new_ulen = row_end-ustart-1;
       new_left =  MAX(-new_ulen, diag);
       new_right = MIN(right_dist, (GtWord)col_end);

       rpoint = evaluateaffineDBcrosspoints(spacemanager, Diagcolumn,
                                        scorehandler, edge, from_edge,Affine_D,
                                        rowoffset,coloffset, useq, ustart,
                                        new_ulen, vseq, vstart, col_end,
                                        new_left, new_right);
       if (col_start + 1 <= vlen)
       {
         Diagcolumn[col_start+1].val_R.last_type = rpoint.edge;
         Diagcolumn[col_start+1].val_D.last_type = rpoint.edge;
         Diagcolumn[col_start+1].val_I.last_type = rpoint.edge;
       }
       if (rpoint.idx == lastrpoint.idx)
         lastrpoint = rpoint;
       break;
     case Affine_I:
       new_ulen = row_end-ustart;
       new_vlen = col_end-1;
       new_left = MAX(left_dist,
                 -(GtWord)new_ulen);

       new_right = MIN(diag,new_vlen);
       rpoint = evaluateaffineDBcrosspoints(spacemanager, Diagcolumn,
                               scorehandler,edge, from_edge, Affine_I,rowoffset,
                               coloffset, useq, ustart, new_ulen, vseq, vstart,
                               new_vlen, new_left, new_right);

       Diagcolumn[col_start].val_I.last_type = rpoint.edge;
       break;
     default:
       gt_assert(false);
     }
   }
   else if (cp_type == Affine_D)
   {
     Diagcolumn[1].val_I.last_type = Affine_R;
     Diagcolumn[1].val_D.last_type = Affine_R;
     Diagcolumn[1].val_R.last_type = Affine_R;
     Diagcolumn[0].val_R.currentrowindex = rowoffset;
     Diagcolumn[0].val_R.last_type = from_edge;
   }
  return lastrpoint;
}

/* calculating alignment in linear space within a specified diagonal band
 * with affine gapcosts */
static void gt_calc_diagonalbandaffinealign(GtLinspaceManagement *spacemanager,
                                            const GtScoreHandler *scorehandler,
                                            GtAlignment *align,
                                            const GtUchar *useq,
                                            GtUword ustart, GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart, GtUword vlen,
                                            GtWord left_dist,
                                            GtWord right_dist)
{
  GtAffinealignDPentry *Atabcolumn;
  GtAffineAlignRtabentry *Rtabcolumn;
  GtAffineDiagAlignentry *Diagcolumn;
  GtAffineAlignRnode lastnode;
  GtUword idx, gap_extension;
  gt_assert(align && scorehandler);

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false); /* no global alignment */
  }

  gt_linspace_management_set_ulen(spacemanager, ulen);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);

  if (ulen == 0UL)
  {
    (void) gt_reconstructalignment_trivial_insertion(align, vlen,
                                                     gap_extension);
    return;

  }
  else if (vlen == 0UL)
  {
    (void) gt_reconstructalignment_trivial_deletion(align, ulen,
                                                    gap_extension);
     return;
  }
  if (gt_linspace_management_checksquare(spacemanager, ulen, vlen,
                                         sizeof (*Atabcolumn),
                                         sizeof (*Rtabcolumn)))
  {
    (void) gt_diagonalbandalign_affinegapcost_in_square_space_generic(
                                                         spacemanager,
                                                         scorehandler, align,
                                                         useq, ustart, ulen,
                                                         vseq, vstart, vlen,
                                                         left_dist, right_dist);
    return;
  }
  gt_linspace_management_check(spacemanager, MIN(right_dist-left_dist,ulen),
                               vlen, sizeof (*Atabcolumn), sizeof (*Rtabcolumn),
                               sizeof (*Diagcolumn));
  Diagcolumn = gt_linspace_management_get_crosspointTabspace(spacemanager);

  /* initialize Diagcolumn */
  for (idx = 0; idx <= vlen; idx++)
  {
    Diagcolumn[idx].val_R =
                      (GtDiagAlignentry) {GT_UWORD_MAX, GT_UWORD_MAX, Affine_X};
    Diagcolumn[idx].val_D =
                      (GtDiagAlignentry) {GT_UWORD_MAX, GT_UWORD_MAX, Affine_X};
    Diagcolumn[idx].val_I =
                      (GtDiagAlignentry) {GT_UWORD_MAX, GT_UWORD_MAX, Affine_X};
  }

  lastnode = evaluateaffineDBcrosspoints(spacemanager, Diagcolumn, scorehandler,
                                         Affine_X, Affine_X, Affine_X, 0, 0,
                                         useq, ustart, ulen, vseq, vstart, vlen,
                                         left_dist, right_dist);
  /* reconstruct alignment */
  gt_reconstructalignment_from_affineDtab(align, Diagcolumn, lastnode.edge,
                                          useq, ulen, vseq, vlen);
}

/* compute alignment with affine gapcosts within a diagonal band */
void gt_diagonalbandalign_affinegapcost_compute_generic(GtLinspaceManagement
                                                        *spacemanager,
                                                        const GtScoreHandler
                                                        *scorehandler,
                                                        GtAlignment *align,
                                                        const GtUchar *useq,
                                                        GtUword ustart,
                                                        GtUword ulen,
                                                        const GtUchar *vseq,
                                                        GtUword vstart,
                                                        GtUword vlen,
                                                        GtWord left_dist,
                                                        GtWord right_dist)
{
  gt_assert(useq  && vseq && spacemanager && scorehandler && align);

  /* set new bounds, if left_dist or right_dist is out of sequence */
  left_dist = MAX(-(GtWord) ulen,left_dist);
  right_dist = MIN((GtWord) vlen,right_dist);
  gt_alignment_set_seqs(align, useq+ustart, ulen, vseq+vstart, vlen);
  gt_calc_diagonalbandaffinealign(spacemanager, scorehandler, align,
                                             useq, ustart, ulen,
                                             vseq, vstart, vlen,
                                             left_dist, right_dist);
}

/* compute alignment with affine gapcosts within a diagonal band */
void gt_diagonalbandalign_affinegapcost_compute(GtLinspaceManagement
                                                *spacemanager,
                                                GtAlignment *align,
                                                const GtUchar *useq,
                                                GtUword ustart, GtUword ulen,
                                                const GtUchar *vseq,
                                                GtUword vstart, GtUword vlen,
                                                GtWord left_dist,
                                                GtWord right_dist,
                                                GtUword matchcost,
                                                GtUword mismatchcost,
                                                GtUword gap_opening,
                                                GtUword gap_extension)
{
  GtScoreHandler *scorehandler = gt_scorehandler_new(matchcost,
                                                     mismatchcost,
                                                     gap_opening,
                                                     gap_extension);

  gt_diagonalbandalign_affinegapcost_compute_generic(spacemanager, scorehandler,
                                                     align, useq, ustart, ulen,
                                                     vseq, vstart, vlen,
                                                     left_dist, right_dist);

  gt_scorehandler_delete(scorehandler);
}

void gt_diagonalbandalign_affinegapcost_check(GT_UNUSED bool forward,
                                              const GtUchar *useq, GtUword ulen,
                                              const GtUchar *vseq, GtUword vlen)
{
  GtUword affine_cost1, affine_cost2, affine_cost3,
          matchcost = 0, mismatchcost = 1,
          gap_opening = 2, gap_extension = 1;
  GtWord left_dist, right_dist;
  GtAlignment *align;
  GtScoreHandler *scorehandler;
  GtLinspaceManagement *spacemanager;

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

  left_dist = -ulen;
  right_dist = vlen;

  scorehandler = gt_scorehandler_new(matchcost, mismatchcost,
                                     gap_opening, gap_extension);
  gt_scorehandler_plain(scorehandler);
  gt_scorehandler_downcase(scorehandler);
  affine_cost1 = gt_diagonalbandalign_affinegapcost_square_space_distance_only(
                                           useq, 0, ulen, vseq, 0, vlen,
                                           left_dist, right_dist, scorehandler);

  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);
  spacemanager = gt_linspace_management_new();

  gt_calc_diagonalbandaffinealign(spacemanager, scorehandler, align,
                                  useq, 0, ulen, vseq, 0, vlen,
                                  left_dist, right_dist);

  gt_linspace_management_delete(spacemanager);
  affine_cost2 = gt_alignment_eval_with_affine_score(align, true, matchcost,
                                                     mismatchcost,
                                                     gap_opening,
                                                     gap_extension);
  if (affine_cost1 != affine_cost2)
  {
    fprintf(stderr,"gt_diagonalband_affinegapcost_square_space_distance_only = "
            GT_WU" != "GT_WU" = gt_alignment_eval_generic_with_affine_score\n",
            affine_cost1, affine_cost2);

    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  affine_cost3 = diagonalband_linear_affine(useq, 0, ulen, vseq, 0, vlen,
                                           left_dist, right_dist, scorehandler);
  if (affine_cost3 != affine_cost2)
  {
    fprintf(stderr,"diagonalband_linear_affine = "GT_WU
            " != "GT_WU" = gt_alignment_eval_generic_with_affine_score\n",
            affine_cost3, affine_cost2);

    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_scorehandler_delete(scorehandler);
  gt_alignment_delete(align);
}
