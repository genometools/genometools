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
#include "core/assert_api.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/error.h"
#include "core/types_api.h"
#include "core/divmodmul.h"
#include "core/unused_api.h"
#include "extended/diagonalbandalign.h"
#include "extended/linspace_management.h"
#include "extended/reconstructalignment.h"
#include "match/squarededist.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

static void diagonalband_fillDPtab_in_square_space(GtUword **EDtabcolumn,
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
  GtUword i,j, val, low_row, high_row, gapcost;
  gt_assert(EDtabcolumn && scorehandler);

  low_row = 0;
  high_row = -left_dist;
  gapcost = gt_scorehandler_get_gapscore(scorehandler);

  /* first column */
  EDtabcolumn[0][0] = 0;

  for (i = 1; i <= high_row; i++)
  {
      EDtabcolumn[i][0] = add_safe_umax(EDtabcolumn[i-1][0], gapcost);
  }
  for (; i <= ulen; i++)
  {
      EDtabcolumn[i][0] = GT_UWORD_MAX; /* invalid value */
  }

  /* next columns */
  for (j = 1; j <= vlen; j++)
  {
    /* below diagonal band*/
    for (i = 0; i <= low_row; i++)
    {
      if (j <= right_dist)
      {
        EDtabcolumn[i][j] = add_safe_umax(EDtabcolumn[i][j-1], gapcost);
      }
      else{
        EDtabcolumn[i][j] = GT_UWORD_MAX;
      }
    }
    if (j > right_dist)
      low_row ++;
    if (high_row < ulen)
      high_row ++;

    /* diagonaldband */
    for (; i <= high_row; i++)
    {
      EDtabcolumn[i][j] = add_safe_umax(EDtabcolumn[i][j-1], gapcost);

      val = add_safe_umax(EDtabcolumn[i-1][j-1],
                          gt_scorehandler_get_replacement(scorehandler,
                                    useq[ustart+i-1],vseq[vstart+j-1]));
      if (val <= EDtabcolumn[i][j])
      {
        EDtabcolumn[i][j] = val;
      }

      if ((val = add_safe_umax(EDtabcolumn[i-1][j],gapcost))
                                        < EDtabcolumn[i][j])
      {
        EDtabcolumn[i][j] = val;
      }
    }
    /* above diagonal band */
    for (; i <= ulen; i++)
      EDtabcolumn[i][j] = GT_UWORD_MAX;
  }
}
/* calculate only distance with diagonalband in square space */
static GtUword diagonalband_squarespace_distance_only(const GtUchar *useq,
                                                      GtUword ustart,
                                                      GtUword ulen,
                                                      const GtUchar *vseq,
                                                      GtUword vstart,
                                                      GtUword vlen,
                                                      GtWord left_dist,
                                                      GtWord right_dist,
                                                      GtUword matchcost,
                                                      GtUword mismatchcost,
                                                      GtUword gapcost)
{
  GtUword **EDtabcolumn, distance = GT_UWORD_MAX;
  GtScoreHandler *scorehandler;

   if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }
  gt_array2dim_malloc(EDtabcolumn, (ulen+1), (vlen+1));
  scorehandler = gt_scorehandler_new(matchcost, mismatchcost, 0, gapcost);
  gt_scorehandler_plain(scorehandler);
  gt_scorehandler_downcase(scorehandler);
  diagonalband_fillDPtab_in_square_space(EDtabcolumn, useq, ustart, ulen,
                                         vseq, vstart, vlen,
                                         left_dist, right_dist, scorehandler);

  distance = EDtabcolumn[ulen][vlen];
  gt_array2dim_delete(EDtabcolumn);
  gt_scorehandler_delete(scorehandler);
  return distance;
}

/* creating alignment with diagonalband in square space */
GtUword gt_diagonalbandalignment_in_square_space_generic(GtLinspaceManagement
                                                         *space,
                                                         GtAlignment *align,
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
  GtUword **EDtabcolumn, distance;

  gt_assert(align && scorehandler);

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }

  if (space == NULL)
  {
    /*use it in normally case*/
    gt_array2dim_malloc(EDtabcolumn, (ulen+1), (vlen+1));
  }
  else
  {
    /*use it in lineraspace context*/
    EDtabcolumn = gt_linspace_management_change_to_square(space, ulen, vlen);
  }

  diagonalband_fillDPtab_in_square_space(EDtabcolumn, useq, ustart, ulen,
                                         vseq, vstart, vlen,
                                         left_dist, right_dist, scorehandler);

  distance = EDtabcolumn[ulen][vlen];
  /* reconstruct alignment from 2dimarray EDtabcolumn */
  gt_reconstructalignment_from_EDtab(align, EDtabcolumn, useq, ustart, ulen,
                                     vseq, vstart, vlen, scorehandler);

  if (space == NULL)
  {
    gt_array2dim_delete(EDtabcolumn);
  }
  return distance;
}

/* creating alignment with diagonalband in square space
   with constant cost values */
GtUword gt_diagonalbandalignment_in_square_space(GtLinspaceManagement
                                                 *spacemanager,
                                                 GtAlignment *alignment,
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
                                                 GtUword gapcost)
{
  GtUword distance;
  GtScoreHandler *scorehandler;

  scorehandler = gt_scorehandler_new(matchcost, mismatchcost, 0, gapcost);
  distance = gt_diagonalbandalignment_in_square_space_generic(spacemanager,
                                                           alignment,
                                                           useq, ustart, ulen,
                                                           vseq, vstart, vlen,
                                                           left_dist,
                                                           right_dist,
                                                           scorehandler);
  gt_scorehandler_delete(scorehandler);
  return distance;
}

static void evaluate_DBcrosspoints_from_2dimtab(GtUword **E,
                                                GtDiagAlignentry *Dtab,
                                                const GtUchar *useq,
                                                GtUword ustart,
                                                GtUword ulen,
                                                const GtUchar *vseq,
                                                GtUword vstart,
                                                GtUword vlen,
                                                GtUword rowoffset,
                                                const GtScoreHandler
                                                *scorehandler)
{
  GtUword idx, jdx, gapcost;

  gt_assert(E);
  gapcost = gt_scorehandler_get_gapscore(scorehandler);

  idx = ulen;
  jdx = vlen;

  while (jdx > 0 || idx > 0)
  {
    if (idx > 0 && jdx > 0 && E[idx][jdx] == E[idx-1][jdx-1] +
                                  gt_scorehandler_get_replacement(scorehandler,
                                       useq[ustart+idx-1], vseq[vstart+jdx-1]))
    {
      if (jdx == vlen)
        Dtab[vlen].currentrowindex = idx + rowoffset;

      Dtab[jdx].last_type = Linear_R;
      idx--;
      jdx--;
      Dtab[jdx].currentrowindex = idx + rowoffset;
    }
    else if (idx > 0 && E[idx][jdx] == E[idx-1][jdx] + gapcost)
    {
      if (jdx == vlen)
        Dtab[vlen].currentrowindex = idx + rowoffset;

      Dtab[jdx].last_type = Linear_D;
      idx--;
      Dtab[jdx].currentrowindex = idx + rowoffset;
    }
    else if (jdx > 0 && E[idx][jdx] == E[idx][jdx-1] + gapcost)
    {
      if (jdx == vlen)
        Dtab[vlen].currentrowindex = idx + rowoffset;

      Dtab[jdx].last_type = Linear_I;
      jdx--;
      Dtab[jdx].currentrowindex = idx + rowoffset;
    }
    else
      gt_assert(false);
  }

}

/*create DBcrosspointtab to combine square calculating with linear calculating*/
GT_UNUSED  static void dtab_in_square_space(GtLinspaceManagement *spacemanager,
                                            GtDiagAlignentry *Dtab,
                                            const GtUchar *useq,
                                            GtUword ustart,
                                            GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart,
                                            GtUword vlen,
                                            GtWord left_dist,
                                            GtWord right_dist,
                                            GtUword rowoffset,
                                            LinearAlignEdge edge,
                                            const GtScoreHandler *scorehandler)
{
  GtUword **EDtabcolumn;
  gt_assert(Dtab && scorehandler && spacemanager);

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }

  EDtabcolumn = gt_linspace_management_change_to_square(spacemanager,ulen,vlen);
  diagonalband_fillDPtab_in_square_space(EDtabcolumn, useq, ustart, ulen,
                                         vseq, vstart, vlen, left_dist,
                                         right_dist, scorehandler);

  evaluate_DBcrosspoints_from_2dimtab(EDtabcolumn, Dtab, useq, ustart, ulen,
                                      vseq, vstart, vlen, rowoffset,
                                      scorehandler);
  Dtab[0].last_type = edge;
}

/* calculate only distance with diagonalband in linear space O(n) */

static GtUword diagonalband_linear_distance_only(const GtUchar *useq,
                                                 GtUword ustart,
                                                 GtUword ulen,
                                                 const GtUchar *vseq,
                                                 GtUword vstart,
                                                 GtUword vlen,
                                                 GtWord left_dist,
                                                 GtWord right_dist,
                                                 GtUword matchcost,
                                                 GtUword mismatchcost,
                                                 GtUword gapcost)
{
  GtUword distance, colindex, rowindex, low_row, high_row, width, val,
        *EDtabcolumn, northwestEDtabentry, westEDtabentry = GT_UWORD_MAX;
  bool last_row = false;

  distance = GT_UWORD_MAX;

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    return GT_UWORD_MAX;
  }

  width = right_dist - left_dist + 1;
  EDtabcolumn = gt_malloc(sizeof(*EDtabcolumn) * width);

  low_row = 0;
  high_row = -left_dist;
  EDtabcolumn[low_row] = 0;

  for (rowindex = low_row+1; rowindex <= high_row; rowindex ++)
  {
    EDtabcolumn[rowindex-low_row] = EDtabcolumn[rowindex-low_row-1] + gapcost;
  }
  if (high_row == ulen)
    last_row = true;
  for (colindex = 1; colindex <= vlen; colindex++)
  {
    northwestEDtabentry = EDtabcolumn[0];

    if (colindex > right_dist)
    {
      if (low_row != high_row)
        westEDtabentry = EDtabcolumn[1];
      low_row++;
    }
    else
      westEDtabentry = EDtabcolumn[0];

   if (high_row < ulen)
      high_row ++;
   if (!last_row && low_row == high_row)
        westEDtabentry = GT_UWORD_MAX;
    EDtabcolumn[0] = add_safe_umax(westEDtabentry, gapcost);

    if (low_row > 0 )
    {
      if ((val = add_safe_umax(northwestEDtabentry,
                              (tolower((int)useq[ustart+low_row-1]) ==
                              tolower((int)vseq[vstart+colindex-1])?
                              matchcost : mismatchcost)))
          < EDtabcolumn[0])
      EDtabcolumn[0] = val;
    }
    for (rowindex = low_row + 1; rowindex <= high_row; rowindex++)
    {
      northwestEDtabentry = westEDtabentry;
      if (!last_row && rowindex == high_row)
        westEDtabentry = GT_UWORD_MAX;
      else if (low_row > 0)
        westEDtabentry = EDtabcolumn[rowindex-low_row+1];
      else
        westEDtabentry = EDtabcolumn[rowindex-low_row];

      if (rowindex == ulen)
        last_row = true;
      EDtabcolumn[rowindex-low_row] = add_safe_umax(westEDtabentry, gapcost);

      val = add_safe_umax(northwestEDtabentry,
                        (tolower((int)useq[ustart+rowindex-1]) ==
                         tolower((int)vseq[vstart+colindex-1]) ?
                         matchcost : mismatchcost));
      if (val <= EDtabcolumn[rowindex-low_row])
        EDtabcolumn[rowindex-low_row] = val;

      if ((val = add_safe_umax(EDtabcolumn[rowindex-low_row-1], gapcost))
                                       <= EDtabcolumn[rowindex-low_row])
        EDtabcolumn[rowindex-low_row] = val;
    }
  }

  distance = EDtabcolumn[high_row-low_row];
  gt_free(EDtabcolumn);

  return distance;
}

static void firstDBtabcolumn(GtUword *EDtabcolumn,
                             GtUword *Rtabcolumn,
                             GtDiagAlignentry *Diagcolumn,
                             LinearAlignEdge edge,
                             const GtWord offset,
                             GtWord left_dist,
                             GtWord right_dist,
                             GtUword gapcost)
{
  GtUword rowindex, low_row, high_row;
  GtWord diag;

  diag = GT_DIV2(left_dist + right_dist);
  low_row = 0;
  high_row = -left_dist;

  EDtabcolumn[low_row] = 0;
  if (diag == 0)
  {
    Diagcolumn[0].last_type = edge;
    Diagcolumn[0].lastcpoint = GT_UWORD_MAX;
    Diagcolumn[0].currentrowindex = 0 + offset;
    Rtabcolumn[0] = 0;
  }
  else
  {
    Rtabcolumn[low_row] = GT_UWORD_MAX;
  }

  for (rowindex = low_row+1; rowindex <= high_row; rowindex++)
  {
    EDtabcolumn[rowindex-low_row] = EDtabcolumn[rowindex-low_row-1] + gapcost;

    if (diag == -(GtWord)rowindex)
    {
      Diagcolumn[0].last_type = Linear_D;
      Diagcolumn[0].lastcpoint = GT_UWORD_MAX;
      Diagcolumn[0].currentrowindex = rowindex + offset;
      Rtabcolumn[rowindex-low_row] = 0;
    }
    else
    {
      Rtabcolumn[rowindex-low_row] = Rtabcolumn[rowindex-low_row-1];
    }
  }
}

static inline void set_linear_DiagentryRtabentry(LinearAlignEdge edge,
                                              GtWord diag, GtUword colindex,
                                              GtUword rowindex,
                                              GtUword offset,
                                              GtDiagAlignentry *Diagcolumnentry,
                                              GtUword *Rtabcolumnentry,
                                              GtUword Rtabentry_from)
{
    if (diag == (GtWord)colindex - (GtWord)rowindex)
    {
      Diagcolumnentry->last_type = edge;
      Diagcolumnentry->lastcpoint = Rtabentry_from;
      Diagcolumnentry->currentrowindex = rowindex + offset;
      *Rtabcolumnentry = colindex;
    }
    else
    {
      *Rtabcolumnentry = Rtabentry_from;
    }
}

/* calculate all E- and Rtabcolumns, store crosspoints in  Diagcolumn,
   return lastcrosspoint from optimal path */
static GtUword evaluateallDBtabcolumns(GtLinspaceManagement *spacemanager,
                                       GtDiagAlignentry *Diagcolumn,
                                       const GtScoreHandler *scorehandler,
                                       LinearAlignEdge edge,
                                       GtWord offset,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen,
                                       GtWord left_dist,
                                       GtWord right_dist)
{
  GtUword gapcost, colindex, rowindex, val, *EDtabcolumn, *Rtabcolumn,
          northwestEDtabentry, westEDtabentry = GT_UWORD_MAX,
          northwestRtabentry, westRtabentry = GT_UWORD_MAX,
          low_row, high_row; /*lowest and highest row between a diagonal band*/
  GtWord diag;
  bool last_row = false;

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }

  EDtabcolumn = gt_linspace_management_get_valueTabspace(spacemanager);
  Rtabcolumn = gt_linspace_management_get_rTabspace(spacemanager);

  gapcost = gt_scorehandler_get_gapscore(scorehandler);
  diag = GT_DIV2(left_dist + right_dist);
  low_row = 0;
  high_row = -left_dist;

 /* first column */
  firstDBtabcolumn(EDtabcolumn, Rtabcolumn, Diagcolumn, edge, offset,
                   left_dist, right_dist, gapcost);
   if (high_row == ulen)
     last_row = true;
  /* next columns */
  for (colindex = 1; colindex <= vlen; colindex++)
  {
    northwestEDtabentry = EDtabcolumn[0];
    northwestRtabentry = Rtabcolumn[0];

    if (colindex > right_dist)
    {
      if (low_row != high_row)
      {
        westEDtabentry = EDtabcolumn[1];
        westRtabentry = Rtabcolumn[1];
      }
      low_row++;
    }
    else
    {
      westEDtabentry = EDtabcolumn[0];
      westRtabentry = Rtabcolumn[0];
    }
    if (high_row < ulen)
      high_row ++;
    if (!last_row && low_row == high_row)
    {/* prev is outside of diagonalband*/
      westEDtabentry = GT_UWORD_MAX;
      westRtabentry = GT_UWORD_MAX;
    }

    EDtabcolumn[0] = add_safe_umax(westEDtabentry, gapcost);
    edge = Linear_I;

    /* replacement possible for 0-entry */
    if (low_row > 0 )
    {
      val = add_safe_umax(northwestEDtabentry,
                         gt_scorehandler_get_replacement(scorehandler,
                         useq[ustart+low_row-1], vseq[vstart+colindex-1]));
      if (val <= EDtabcolumn[0])
      {
        edge = Linear_R;
        EDtabcolumn[0] = val;
      }
    }

    if (edge == Linear_R)
    {
      set_linear_DiagentryRtabentry(edge, diag, colindex, low_row, offset,
                                    &Diagcolumn[colindex],&Rtabcolumn[0],
                                    northwestRtabentry);
    } else
    {
      gt_assert(edge == Linear_I);
      set_linear_DiagentryRtabentry(edge, diag, colindex, low_row, offset,
                                    &Diagcolumn[colindex],&Rtabcolumn[0],
                                    westRtabentry);
    }

    for (rowindex = low_row + 1; rowindex <= high_row; rowindex++)
    {
      GtUword Rtabentry_from = 0;
      northwestEDtabentry = westEDtabentry;
      northwestRtabentry = westRtabentry;

      if (!last_row && rowindex == high_row)
      {/* prev is outside of diagonalband*/
        westEDtabentry = GT_UWORD_MAX;
        westRtabentry = GT_UWORD_MAX;
      }
      else if (low_row > 0)
      {/* shifted diagonalband*/
        westEDtabentry = EDtabcolumn[rowindex-low_row+1];
        westRtabentry = Rtabcolumn[rowindex-low_row+1];
      }
      else
      {/* normaly prev*/
        westEDtabentry = EDtabcolumn[rowindex-low_row];
        westRtabentry = Rtabcolumn[rowindex-low_row];
      }

      if (rowindex == ulen)
        last_row = true;
      /* insertion */
      EDtabcolumn[rowindex-low_row] = add_safe_umax(westEDtabentry, gapcost);
      edge = Linear_I;

      /* replacement */
      val = add_safe_umax(northwestEDtabentry,
                          gt_scorehandler_get_replacement(scorehandler,
                          useq[ustart+rowindex-1], vseq[vstart+colindex-1]));

      if (val <= EDtabcolumn[rowindex-low_row])
      {
        EDtabcolumn[rowindex-low_row] = val;
        edge = Linear_R;
      }
      /* deletion */
      if ((val = add_safe_umax(EDtabcolumn[rowindex-low_row-1], gapcost))
                                       < EDtabcolumn[rowindex-low_row])
      {
        EDtabcolumn[rowindex-low_row] = val;
        edge = Linear_D;
      }

      if (edge == Linear_R)
      {
        Rtabentry_from = northwestRtabentry;
      } else
      {
        if (edge == Linear_D)
        {
          Rtabentry_from = Rtabcolumn[rowindex-low_row-1];
        } else
        {
          gt_assert(edge == Linear_I);
          Rtabentry_from = westRtabentry;
        }
      }
      set_linear_DiagentryRtabentry(edge, diag, colindex, rowindex, offset,
                           &Diagcolumn[colindex],&Rtabcolumn[rowindex-low_row],
                           Rtabentry_from);
    }
  }
  return Rtabcolumn[high_row-low_row];
}

/* calculate crosspoint realting to diagonal in recursive way */
static void evaluateDBcrosspoints(GtLinspaceManagement *spacemanager,
                                  GtDiagAlignentry *Diagcolumn,
                                  const GtScoreHandler *scorehandler,
                                  LinearAlignEdge edge,
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
  GtUword idx, prevcpoint, cpoint, ctemp, new_ulen;
  GtWord new_left, new_right, diag = GT_DIV2(left_dist+right_dist);
  GtDiagAlignentry dtemp;

  if (ulen == 0)
  {
    for (idx = 1; idx <=vlen; idx++)
    {
      Diagcolumn[idx].currentrowindex = rowoffset;
      Diagcolumn[idx].last_type = Linear_I;
    }
    Diagcolumn[0].currentrowindex = rowoffset;
    Diagcolumn[0].last_type = edge;
    return;
  }

  if (vlen == 0)
  {
    Diagcolumn[0] = (GtDiagAlignentry) {ulen, edge};
    return;
  }

  if (gt_linspace_management_checksquare(spacemanager,ulen,vlen,
                                         sizeof (GtUword),sizeof(GtUword)))
  {
    dtab_in_square_space(spacemanager, Diagcolumn,
                         useq, ustart, ulen,
                         vseq, vstart, vlen,
                         left_dist, right_dist,
                         rowoffset, edge, scorehandler);
    return;
  }

  cpoint = evaluateallDBtabcolumns (spacemanager, Diagcolumn, scorehandler,edge,
                                    rowoffset, useq, ustart, ulen,
                                    vseq, vstart, vlen,
                                    left_dist, right_dist);

  /* if no crosspoint is found */
  if (cpoint == GT_UWORD_MAX)
  {
    if (diag < 0)
    {
      return evaluateDBcrosspoints(spacemanager,Diagcolumn, scorehandler, edge,
                                   rowoffset, coloffset, useq, ustart, ulen,
                                   vseq, vstart, vlen,
                                   diag+1, right_dist);
    }

    else if (diag > 0)
      return evaluateDBcrosspoints(spacemanager, Diagcolumn, scorehandler, edge,
                                   rowoffset, coloffset, useq, ustart, ulen,
                                   vseq, vstart, vlen,
                                   left_dist, diag-1);
    else
    {
      gt_assert(false); /* there have to be an crosspoint */
    }
  }

  /* exception, if last crosspoint != (m+1)entry, bottom right corner */
  if (cpoint != vlen)
  {
    if (diag + ((GtWord)ulen-(GtWord)vlen) > 0)
    {
      dtemp = Diagcolumn[cpoint];
      new_left = MAX((GtWord)left_dist-diag+1,
                    -((GtWord)ulen-((GtWord)Diagcolumn[cpoint].currentrowindex+1
                    -(GtWord)rowoffset)));
      new_right = 0;
      new_ulen =  ulen - (Diagcolumn[cpoint].currentrowindex+1-rowoffset);
      evaluateDBcrosspoints(spacemanager, Diagcolumn+cpoint, scorehandler,
                            Linear_D, Diagcolumn[cpoint].currentrowindex+1,
                            coloffset+cpoint, useq,
                            Diagcolumn[cpoint].currentrowindex+1, new_ulen,
                            vseq, vstart+cpoint, vlen-cpoint,
                            new_left, new_right);

    Diagcolumn[cpoint] = dtemp;
    }
    else
    {
      new_left = -1;
      new_right =  MIN((GtWord)right_dist-((GtWord)diag)-1,
                      ((GtWord)vlen-(GtWord)cpoint-1));
      new_ulen = ulen - (Diagcolumn[cpoint].currentrowindex-rowoffset);
      evaluateDBcrosspoints(spacemanager,Diagcolumn+cpoint+1,scorehandler,
                            Linear_I, Diagcolumn[cpoint].currentrowindex,
                            coloffset+cpoint+1, useq,
                            Diagcolumn[cpoint].currentrowindex, new_ulen,
                            vseq, vstart+cpoint+1, vlen-cpoint-1,
                            new_left, new_right);
    }
  }

  /* look at all 'normally' crosspoints,
     call segments between pairwise crosspoints*/
  while (Diagcolumn[cpoint].lastcpoint != GT_UWORD_MAX)
  {
    prevcpoint = cpoint;
    if (prevcpoint == 0)
      break;

    cpoint = Diagcolumn[cpoint].lastcpoint;
    ctemp =  Diagcolumn[cpoint].lastcpoint;
    if (Diagcolumn[prevcpoint].last_type == Linear_R ||
       ((Diagcolumn[prevcpoint].last_type == Linear_I) &&
       (prevcpoint-cpoint == 1)))
    {
      continue;
    }
    else if (Diagcolumn[prevcpoint].last_type == Linear_D)
    {
      new_left = -1;
      new_right = MIN(right_dist-((GtWord)diag)-1,
                     (GtWord)prevcpoint-(GtWord)cpoint-1);
      new_ulen = Diagcolumn[prevcpoint].currentrowindex-
                 Diagcolumn[cpoint].currentrowindex-1;

      evaluateDBcrosspoints(spacemanager,Diagcolumn+cpoint+1,scorehandler,
                   Linear_I, Diagcolumn[cpoint].currentrowindex,
                   coloffset+cpoint+1, useq, Diagcolumn[cpoint].currentrowindex,
                   new_ulen, vseq, vstart + cpoint+1, prevcpoint-cpoint-1,
                   new_left, new_right);
    }
    else if (Diagcolumn[prevcpoint].last_type == Linear_I)
    {
      dtemp = Diagcolumn[cpoint];
      new_ulen = Diagcolumn[prevcpoint].currentrowindex-
                 Diagcolumn[cpoint].currentrowindex-1;
      new_left = MAX(left_dist-diag+1,-new_ulen);
      new_right = 0;

      evaluateDBcrosspoints(spacemanager,Diagcolumn+cpoint,scorehandler,
                            Linear_D, Diagcolumn[cpoint].currentrowindex+1,
                            coloffset+cpoint, useq,
                            Diagcolumn[cpoint].currentrowindex+1, new_ulen,
                            vseq, vstart + cpoint, prevcpoint-1-cpoint,
                            new_left, new_right);
      Diagcolumn[cpoint] = dtemp;
    }
    else
    {
      /* if (Diagcolumn[cpoint].last_type == Linear_X), never reach this line */
      gt_assert(false);
    }
    Diagcolumn[cpoint].lastcpoint  = ctemp;
  }

  /* exception, if first crosspoint != 0-entry, upper left corner */
  if (vstart-coloffset != cpoint)
  {
    if (Diagcolumn[cpoint].last_type == Linear_D)
    {
      new_ulen = Diagcolumn[cpoint].currentrowindex-ustart-1;
      new_left =  MAX(diag, -new_ulen);
      new_right = MIN(right_dist, (GtWord)cpoint);

      evaluateDBcrosspoints(spacemanager, Diagcolumn, scorehandler,
                            edge, rowoffset, coloffset,
                            useq, ustart, new_ulen,
                            vseq, vstart, cpoint,
                            new_left, new_right);

    }
    else if (Diagcolumn[cpoint].last_type == Linear_I)
    {
      new_left = MAX(left_dist,
                 -((GtWord)Diagcolumn[cpoint].currentrowindex-(GtWord)ustart));
      new_right = MIN((GtWord)cpoint-1, diag);
      evaluateDBcrosspoints(spacemanager, Diagcolumn, scorehandler,
                            edge, rowoffset, coloffset, useq, ustart,
                            Diagcolumn[cpoint].currentrowindex-ustart,
                            vseq, vstart, cpoint-1,
                            new_left, new_right);
    }
    else if (Diagcolumn[cpoint].last_type == Linear_R ||
             Diagcolumn[cpoint].last_type  == Linear_X)
    {
      gt_assert(false); /* this cant be the first crosspoint or
                           it have to be 0-entry */
    }
  }
}

/* calculating alignment in linear space within a specified diagonal band */
static void gt_calc_diagonalbandalign(GtLinspaceManagement *spacemanager,
                                      const GtScoreHandler *scorehandler,
                                      GtAlignment *align,
                                      const GtUchar *useq,
                                      GtUword ustart, GtUword ulen,
                                      const GtUchar *vseq,
                                      GtUword vstart, GtUword vlen,
                                      GtWord left_dist,
                                      GtWord right_dist)
{
  GtDiagAlignentry *Diagcolumn;
  GtUword idx, gapcost, *EDtabcolumn, *Rtabcolumn;

  gt_assert(align && spacemanager && scorehandler);

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false); /* no global alignment */
  }

  gt_linspace_management_set_ulen(spacemanager,ulen);
  gapcost = gt_scorehandler_get_gapscore(scorehandler);
  if (ulen == 0UL)
  {
    (void) gt_reconstructalignment_trivial_insertion(align, vlen, gapcost);
    return;
  }
  else if (vlen == 0UL)
  {
    (void) gt_reconstructalignment_trivial_deletion(align, ulen, gapcost);
    return;
  }
  else if (gt_linspace_management_checksquare(spacemanager, ulen, vlen,
                                              sizeof (*EDtabcolumn),
                                              sizeof (*Rtabcolumn)))
  {
    (void) gt_diagonalbandalignment_in_square_space_generic(spacemanager,align,
                                                            useq, ustart, ulen,
                                                            vseq, vstart, vlen,
                                                            left_dist,
                                                            right_dist,
                                                            scorehandler);
    return;
  }

  gt_linspace_management_check(spacemanager, MIN(right_dist-left_dist,ulen),
                               vlen,
                               sizeof (*EDtabcolumn),
                               sizeof (*Rtabcolumn),
                               sizeof (*Diagcolumn));

  Diagcolumn = gt_linspace_management_get_crosspointTabspace(spacemanager);
  /* initialize Diagcolumn */
  for (idx = 0; idx <= vlen; idx++)
  {
    Diagcolumn[idx].lastcpoint = GT_UWORD_MAX;
    Diagcolumn[idx].currentrowindex = GT_UWORD_MAX;
  }

  evaluateDBcrosspoints(spacemanager, Diagcolumn, scorehandler,
                        Linear_X, 0, 0,  useq, ustart, ulen,
                        vseq, vstart, vlen, left_dist, right_dist);

  gt_reconstructalignment_from_Dtab(align, Diagcolumn, ulen, vlen);
}

/* compute alignment within a diagonal band */
void gt_diagonalbandalign_compute_generic(GtLinspaceManagement *spacemanager,
                                          const GtScoreHandler *scorehandler,
                                          GtAlignment *align,
                                          const GtUchar *useq,
                                          GtUword ustart, GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vstart, GtUword vlen,
                                          GtWord left_dist,
                                          GtWord right_dist)
{

  gt_assert(useq && vseq && align && spacemanager && scorehandler);

  /* set new bounds, if left_dist or right_dist is out of sequence */
  left_dist = MAX(-(GtWord) ulen, left_dist);
  right_dist = MIN((GtWord) vlen, right_dist);

  gt_alignment_set_seqs(align,useq+ustart, ulen, vseq+vstart, vlen);
  gt_calc_diagonalbandalign(spacemanager, scorehandler, align,
                            useq, ustart, ulen, vseq, vstart, vlen,
                            left_dist, right_dist);
}

void gt_diagonalbandalign_compute(GtLinspaceManagement *spacemanager,
                                  GtAlignment *align,
                                  const GtUchar *useq,
                                  GtUword ustart, GtUword ulen,
                                  const GtUchar *vseq,
                                  GtUword vstart, GtUword vlen,
                                  GtWord left_dist,
                                  GtWord right_dist,
                                  GtUword matchcost,
                                  GtUword mismatchcost,
                                  GtUword gapcost)
{
  GtScoreHandler *scorehandler;
  gt_assert(spacemanager && align);

  scorehandler = gt_scorehandler_new(matchcost, mismatchcost, 0, gapcost);

  gt_diagonalbandalign_compute_generic(spacemanager, scorehandler, align,
                                       useq, ustart, ulen,
                                       vseq, vstart, vlen,
                                       left_dist, right_dist);
  gt_scorehandler_delete(scorehandler);
}

void gt_diagonalbandalign_check(GT_UNUSED bool forward,
                                const GtUchar *useq,
                                GtUword ulen,
                                const GtUchar *vseq,
                                GtUword vlen)
{
  GtUword edist1, edist2, edist3, matchcost = 0, mismatchcost = 1, gapcost = 1;
  GtWord left_dist, right_dist;
  GtAlignment *align;
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

  /* set valid left_dist and right_dist for this test*/
  if ((GtWord)vlen-(GtWord)ulen > 0)
    left_dist = 0;
  else
    left_dist = (GtWord)vlen-(GtWord)ulen-1;
  if ((GtWord)vlen-(GtWord)ulen > 0)
    right_dist = (GtWord)vlen-(GtWord)ulen+2;
  else
    right_dist = 0;

  scorehandler = gt_scorehandler_new(matchcost, mismatchcost, 0, gapcost);
  gt_scorehandler_plain(scorehandler);
  gt_scorehandler_downcase(scorehandler);
  edist1 = diagonalband_linear_distance_only(useq, 0, ulen, vseq, 0, vlen,
                                             left_dist, right_dist, matchcost,
                                             mismatchcost, gapcost);

  edist2 = diagonalband_squarespace_distance_only(useq, 0, ulen, vseq, 0, vlen,
                                                  left_dist, right_dist,
                                                  matchcost, mismatchcost,
                                                  gapcost);

  if (edist1 != edist2)
  {
    fprintf(stderr,"diagonalband_linear_distance_only = "GT_WU" != "GT_WU
              " = diagonalband_squarespace_distance_only\n", edist1, edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  spacemanager = gt_linspace_management_new();
  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);

  gt_calc_diagonalbandalign(spacemanager, scorehandler, align,
                            useq, 0, ulen, vseq, 0, vlen,
                            left_dist, right_dist);

  edist3 = gt_alignment_eval_with_score(align, true, matchcost,
                                        mismatchcost, gapcost);
  if (edist2 != edist3)
  {
    fprintf(stderr,"diagonalband_squarespace_distance_only = "GT_WU" != "GT_WU
              " = gt_alignment_eval_with_score\n", edist2, edist3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_alignment_delete(align);
  gt_linspace_management_delete(spacemanager);
  gt_scorehandler_delete(scorehandler);
}
