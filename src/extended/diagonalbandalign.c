/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
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
#include "extended/linearalign_utilities.h"
#include "extended/reconstructalignment.h"
#include "match/squarededist.h"

#include "extended/diagonalbandalign.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)
typedef enum {
  Linear_R,
  Linear_D,
  Linear_I,
  Linear_X /* unknown */
} LinearAlignEdge;

typedef struct {
  GtUword lastcpoint, currentrowindex;
  LinearAlignEdge edge;
} Diagentry;

/*reconstruct alignment from crosspoints, crosspoints relating to diagonalband*/
void reconstructalignment_from_Dtab(GtAlignment *align,
                                    const Diagentry *Dtab,
                                    GtUword ulen, GtUword vlen)
{
  GtUword i,j;

  gt_assert(align != NULL && Dtab != NULL);

  for (j = ulen; j > Dtab[vlen].currentrowindex; j--)
  {
    gt_alignment_add_deletion(align);
  }
  for (i = vlen; i > 0; i--) {
    gt_assert(Dtab[i].currentrowindex != GT_UWORD_MAX);
    if (Dtab[i].currentrowindex == Dtab[i-1].currentrowindex + 1)
    {
      if (Dtab[i].edge == Linear_R)
       gt_alignment_add_replacement(align);

      else if (Dtab[i].edge == Linear_D)
      {
         gt_alignment_add_deletion(align);
         gt_alignment_add_insertion(align);
      }
      else if (Dtab[i].edge == Linear_I)
      {
         gt_alignment_add_insertion(align);
         gt_alignment_add_deletion(align);
      }
    }
    else if (Dtab[i].currentrowindex == Dtab[i-1].currentrowindex)
      gt_alignment_add_insertion(align);

    else if (Dtab[i].currentrowindex > Dtab[i-1].currentrowindex)
    {
      if (Dtab[i].edge == Linear_R)
      {
        gt_alignment_add_replacement(align);
        for (j = 0; j < (Dtab[i].currentrowindex -
                         Dtab[i-1].currentrowindex)-1; j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else if (Dtab[i].edge == Linear_I)
      {
        gt_alignment_add_insertion(align);
        for (j = 0; j < (Dtab[i].currentrowindex -
                         Dtab[i-1].currentrowindex); j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else
      {
        for (j = 0; j < (Dtab[i].currentrowindex -
                         Dtab[i-1].currentrowindex)-1; j++)
        {
          gt_alignment_add_deletion(align);
        }
       gt_alignment_add_replacement(align);
      }
    }
  }
  for (j = Dtab[0].currentrowindex; j > 0; j--)
  {
    gt_alignment_add_deletion(align);
  }
}

static void diagonalband_fillDPtab_in_square_space(GtUword **E,
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
  GtUword i,j, val, low_row, high_row;

  gt_assert(E != NULL);
  low_row = 0;
  high_row = -left_dist;

  /* first column */
  E[0][0] = 0;
  for (i = 1; i <= high_row; i++)
  {
      E[i][0] = add_safe_umax(E[i-1][0], gapcost);
  }
  for (; i <= ulen; i++)
  {
      E[i][0] = GT_UWORD_MAX; /* invalid value*/
  }

  /* next columns */
  for (j = 1; j <= vlen; j++)
  {
    /* below diagonal band*/
    for (i = 0; i <= low_row; i++)
    {
      if (j <= right_dist)
      {
        E[i][j] = add_safe_umax(E[i][j-1], gapcost);
      }
      else{
        E[i][j] = GT_UWORD_MAX;
      }
    }
    if ( j > right_dist)
      low_row++;
    if (high_row < ulen)
      high_row ++;

    /* diagonaldband */
    for (; i <= high_row; i++)
    {
      E[i][j] = add_safe_umax(E[i][j-1], gapcost);

      if ((val = add_safe_umax(E[i-1][j-1],(tolower((int)useq[ustart+i-1]) ==
                                           tolower((int)vseq[vstart+j-1]) ?
                                matchcost : mismatchcost)))
          <= E[i][j])
      {
        E[i][j] = val;
      }

      if ((val = add_safe_umax(E[i-1][j],gapcost)) < E[i][j])
      {
        E[i][j] = val;
      }
    }
    /* above diagonal band */
    for (; i <= ulen; i++)
      E[i][j] = GT_UWORD_MAX;
  }
}
/* calculate only distance with diagonalband in square space O(n²) */
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
  GtUword **E, distance = GT_UWORD_MAX;

   if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    return GT_UWORD_MAX;
  }

  gt_array2dim_malloc(E, (ulen+1), (vlen+1));
  diagonalband_fillDPtab_in_square_space(E, useq, ustart, ulen, vseq, vstart,
                                         vlen, left_dist, right_dist, matchcost,
                                         mismatchcost, gapcost);

  distance = E[ulen][vlen];
  gt_array2dim_delete(E);
  return distance;
}

static GtUword diagonalbandalignment_in_square_space(GtAlignment *align,
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
  GtUword **E, distance;

  gt_assert(align != NULL);
  gt_array2dim_malloc(E, (ulen+1), (vlen+1));

  diagonalband_fillDPtab_in_square_space(E, useq, ustart, ulen, vseq, vstart,
                                         vlen, left_dist, right_dist, matchcost,
                                         mismatchcost, gapcost);

  distance = E[ulen][vlen];
  /* reconstruct alignment from 2dimarray E */
  reconstructalignment_from_EDtab(align, E, useq, ustart, ulen, vseq, vstart,
                                  vlen, matchcost, mismatchcost, gapcost);

  gt_array2dim_delete(E);
  return distance;
}

void evaluate_DBcrosspoints_from_2dimtab(GtUword **E,
                                       Diagentry *Dtab,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen,
                                       GtUword matchcost,
                                       GtUword mismatchcost,
                                       GtUword gapcost,
                                       GtUword rowoffset)

{
  GtUword idx, jdx;

  gt_assert(E && Dtab && vlen>0);

  idx = ulen;
  jdx = vlen;

  while (jdx > 0 || idx > 0)
  {
    if (idx > 0 && jdx > 0 && E[idx][jdx] == E[idx-1][jdx-1] +
       (tolower((int)useq[ustart+idx-1]) == tolower((int) vseq[vstart+jdx-1]) ?
                                                     matchcost : mismatchcost))
    {
      if (jdx == vlen)
        Dtab[vlen].currentrowindex = idx + rowoffset;

      Dtab[jdx].edge = Linear_R;
      idx--;
      jdx--;
      Dtab[jdx].currentrowindex = idx + rowoffset;
    }
    else if (idx > 0 && E[idx][jdx] == E[idx-1][jdx] + gapcost)
    {
      if (jdx == vlen)
        Dtab[vlen].currentrowindex = idx + rowoffset;

      Dtab[jdx].edge = Linear_D;
      idx--;
      Dtab[jdx].currentrowindex = idx + rowoffset;
    }
    else if (jdx > 0 && E[idx][jdx] == E[idx][jdx-1] + gapcost)
    {
      if (jdx == vlen)
        Dtab[vlen].currentrowindex = idx + rowoffset;

      Dtab[jdx].edge = Linear_I;
      jdx--;
      Dtab[jdx].currentrowindex = idx + rowoffset;
    }
    else
      gt_assert(false);
  }

}

static void dtab_in_square_space(Diagentry *Dtab,
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
                                 GtUword gapcost,
                                 GtUword rowoffset)
{
  GtUword **E;
  gt_assert(Dtab != NULL);

  gt_array2dim_malloc(E, (ulen+1), (vlen+1));
  diagonalband_fillDPtab_in_square_space(E, useq, ustart, ulen, vseq, vstart,
                                         vlen, left_dist, right_dist, matchcost,
                                         mismatchcost, gapcost);

  evaluate_DBcrosspoints_from_2dimtab(E, Dtab, useq, ustart, ulen, vseq, vstart,
                                      vlen, matchcost, mismatchcost, gapcost,
                                      rowoffset);
  gt_array2dim_delete(E);
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
        *EDtabcolumn, northwestEDtabentry, westEDtabentry;
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
  for (colindex = 1; colindex <= vlen; colindex++)
  {
    northwestEDtabentry = EDtabcolumn[0];

    if (colindex > right_dist)
    {
      westEDtabentry = EDtabcolumn[1];
      low_row++;
    }
    else
      westEDtabentry = EDtabcolumn[0];
    if (high_row < ulen)
      high_row ++;
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
                             Diagentry *Diagcolumn,
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
    Diagcolumn[0].edge = edge;
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
      Diagcolumn[0].edge = Linear_D;
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

/* calculate all E- and Rtabcolumns, store crosspoints in  Diagcolumn,
 * return lastcrosspoint from optimal path */
static GtUword evaluateallDBtabcolumns(GtUword *EDtabcolumn,
                                       GtUword *Rtabcolumn,
                                       Diagentry *Diagcolumn,
                                       LinearAlignEdge edge,
                                       GtWord offset,
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
  GtUword colindex, rowindex, val, northwestEDtabentry, westEDtabentry,
          northwestRtabentry, westRtabentry,
          low_row, high_row; /*lowest and highest row between a diagonal band*/
  GtWord diag;
  bool last_row = false;

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }
  diag = GT_DIV2(left_dist + right_dist);
  low_row = 0;
  high_row = -left_dist;

 /* first column */
  firstDBtabcolumn(EDtabcolumn, Rtabcolumn, Diagcolumn, edge, offset,
                   left_dist, right_dist, gapcost);
  /* next columns */
  for (colindex = 1; colindex <= vlen; colindex++)
  {
    northwestEDtabentry = EDtabcolumn[0];
    northwestRtabentry = Rtabcolumn[0];

    if (colindex > right_dist)
    {
      westEDtabentry = EDtabcolumn[1];
      westRtabentry = Rtabcolumn[1];
      low_row++;
    }
    else
    {
      westEDtabentry = EDtabcolumn[0];
      westRtabentry = Rtabcolumn[0];
    }
    if (high_row < ulen)
      high_row ++;

    EDtabcolumn[0] = add_safe_umax(westEDtabentry, gapcost);
    if (diag == (GtWord)colindex - (GtWord)low_row)
    {
      Diagcolumn[colindex].edge = Linear_I;
      Diagcolumn[colindex].lastcpoint = westRtabentry;
      Diagcolumn[colindex].currentrowindex = low_row + offset;
      Rtabcolumn[0] = colindex;
    }
    else
    {
      Rtabcolumn[0] = westRtabentry;
    }

    /* replacement possible for 0-entry */
    if (low_row > 0 )
    {
      if ((val = add_safe_umax(northwestEDtabentry,
                             (tolower((int)useq[ustart+low_row-1]) ==
                              tolower((int)vseq[vstart+colindex-1])?
                              matchcost : mismatchcost)))
          <= EDtabcolumn[0])
      {
        EDtabcolumn[0] = val;

        if (diag == (GtWord)colindex - (GtWord)low_row)
        {
          Diagcolumn[colindex].edge = Linear_R;
          Diagcolumn[colindex].lastcpoint = northwestRtabentry;
          Diagcolumn[colindex].currentrowindex = low_row + offset;
          Rtabcolumn[0] = colindex;
        }
        else
        {
          Rtabcolumn[0] = northwestRtabentry;
        }
      }
    }
    for (rowindex = low_row + 1; rowindex <= high_row; rowindex++)
    {
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
      if (diag == (GtWord)colindex - (GtWord)rowindex)
      {
        Diagcolumn[colindex].edge = Linear_I;
        Diagcolumn[colindex].lastcpoint = westRtabentry;
        Diagcolumn[colindex].currentrowindex = rowindex+offset;
        Rtabcolumn[rowindex-low_row] = colindex;
      }
      else
      {
        Rtabcolumn[rowindex-low_row] = westRtabentry;
      }
      /* replacement */
      val = add_safe_umax(northwestEDtabentry,
                        (tolower((int)useq[ustart+rowindex-1]) ==
                         tolower((int)vseq[vstart+colindex-1]) ?
                         matchcost : mismatchcost));
      if (val <= EDtabcolumn[rowindex-low_row])
      {
        EDtabcolumn[rowindex-low_row] = val;
        if (diag == (GtWord)colindex - (GtWord)rowindex)
        {
          Diagcolumn[colindex].edge = Linear_R;
          Diagcolumn[colindex].lastcpoint = northwestRtabentry;
          Diagcolumn[colindex].currentrowindex = rowindex + offset;
          Rtabcolumn[rowindex-low_row] = colindex;
        }
        else
        {
          Rtabcolumn[rowindex-low_row] = northwestRtabentry;
        }
      }
      /* deletion */
      if ((val = add_safe_umax(EDtabcolumn[rowindex-low_row-1], gapcost))
                                       < EDtabcolumn[rowindex-low_row])
      {
        EDtabcolumn[rowindex-low_row] = val;
        if (diag == (GtWord)colindex - (GtWord)rowindex)
        {
          Diagcolumn[colindex].edge = Linear_D;
          Diagcolumn[colindex].lastcpoint = Rtabcolumn[rowindex-low_row-1];
          Diagcolumn[colindex].currentrowindex = rowindex + offset;
          Rtabcolumn[rowindex-low_row] = colindex;
        }
        else
        {
          Rtabcolumn[rowindex-low_row] = Rtabcolumn[rowindex-low_row-1];
        }
      }
    }
  }
  /* last crosspoint of optimal path */
  return Rtabcolumn[high_row-low_row];
}

/* calculate crosspoint in recursive way */
static void evaluateDBcrosspoints(GtUword *EDtabcolumn,
                                  GtUword *Rtabcolumn,
                                  Diagentry *Diagcolumn,
                                  LinearAlignEdge edge,
                                  GtUword rowoffset,
                                  GtUword coloffset,
                                  const GtUchar *useq,
                                  GtUword ustart,
                                  GtUword ulen,
                                  GtUword original_ulen,
                                  const GtUchar *vseq,
                                  GtUword vstart,
                                  GtUword vlen,
                                  GtUword original_vlen,
                                  GtWord left_dist,
                                  GtWord right_dist,
                                  GtUword matchcost,
                                  GtUword mismatchcost,
                                  GtUword gapcost)
{
  GtUword idx, prevcpoint, cpoint, ctemp, new_ulen;
  GtWord new_left, new_right, diag = GT_DIV2(left_dist+right_dist);
  Diagentry dtemp;
  if (ulen == 0)
  {
    for (idx = 0; idx <=vlen; idx++)
    {
      Diagcolumn[idx].currentrowindex = rowoffset;
    }
    return;
  }

  if (vlen == 0)
  {
    Diagcolumn[0].currentrowindex = ulen;
    return;
  }

  if ((ulen+1)*(vlen+1) <= (original_ulen+1))
  {
    dtab_in_square_space(Diagcolumn, useq, ustart, ulen, vseq, vstart, vlen,
                         left_dist, right_dist, matchcost, mismatchcost,
                         gapcost, rowoffset);
    return;
  }
  cpoint = evaluateallDBtabcolumns (EDtabcolumn, Rtabcolumn, Diagcolumn, edge,
                              rowoffset, useq, ustart, ulen, vseq, vstart, vlen,
                              left_dist, right_dist,
                              matchcost, mismatchcost, gapcost);

  /* no crosspoint for optimal path*/
  if (cpoint == GT_UWORD_MAX)
  {
    if (diag < 0)
      return evaluateDBcrosspoints(EDtabcolumn, Rtabcolumn,Diagcolumn, edge,
                                   rowoffset, coloffset, useq, ustart, ulen,
                                   original_ulen, vseq, vstart, vlen,
                                   original_vlen, left_dist+1, right_dist,
                                   matchcost,mismatchcost, gapcost);
    else if (diag > 0)
      return evaluateDBcrosspoints(EDtabcolumn, Rtabcolumn, Diagcolumn, edge,
                                   rowoffset, coloffset, useq, ustart, ulen,
                                   original_ulen, vseq, vstart, vlen,
                                   original_vlen, left_dist, right_dist-1,
                                   matchcost, mismatchcost, gapcost);
    else
    {
      gt_assert(false); /* there have to be an crosspoint */
    }
  }

  /* exception, if last cpoint != (m+1)entry */
  if (cpoint != vlen)
  {
    if (diag + ((GtWord)ulen-(GtWord)vlen) > 0)
    {
      dtemp = Diagcolumn[cpoint];
      new_left = MAX((GtWord)left_dist,
                    -((GtWord)ulen-((GtWord)Diagcolumn[cpoint].currentrowindex+1
                    -(GtWord)rowoffset)));
      new_right = 0;
      new_ulen =  ulen - (Diagcolumn[cpoint].currentrowindex+1-rowoffset);
      evaluateDBcrosspoints(EDtabcolumn, Rtabcolumn, Diagcolumn+cpoint,Linear_D,
                         Diagcolumn[cpoint].currentrowindex+1, coloffset+cpoint,
                         useq, Diagcolumn[cpoint].currentrowindex+1, new_ulen,
                         original_ulen, vseq, vstart+cpoint, vlen-cpoint,
                         original_vlen, new_left, new_right,
                         matchcost, mismatchcost, gapcost);

    Diagcolumn[cpoint] = dtemp;
    }
    else
    {
      new_left = -1;
      new_right =  MIN((GtWord)right_dist-((GtWord)diag)-1,
                      ((GtWord)vlen-(GtWord)cpoint-1));
      new_ulen = ulen - (Diagcolumn[cpoint].currentrowindex-rowoffset);
      evaluateDBcrosspoints(EDtabcolumn,Rtabcolumn,Diagcolumn+cpoint+1,Linear_I,
                          Diagcolumn[cpoint].currentrowindex,coloffset+cpoint+1,
                          useq, Diagcolumn[cpoint].currentrowindex, new_ulen,
                          original_ulen, vseq, vstart+cpoint+1, vlen-cpoint-1,
                          original_vlen, new_left, new_right,
                          matchcost, mismatchcost, gapcost);
    }
  }

  /* look at all 'normally' crosspoints */
  while (Diagcolumn[cpoint].lastcpoint != GT_UWORD_MAX)
  {
    prevcpoint = cpoint;
    if (prevcpoint == 0)
      break;

    cpoint = Diagcolumn[cpoint].lastcpoint;
    ctemp =  Diagcolumn[cpoint].lastcpoint;
    if (Diagcolumn[prevcpoint].edge == Linear_R || prevcpoint-cpoint == 1)
    {
      continue;
    }
    else if (Diagcolumn[prevcpoint].edge == Linear_D)
    {
      new_left = -1;
      new_right = MIN(right_dist-((GtWord)diag),
                     (GtWord)prevcpoint-(GtWord)cpoint-1);
      new_ulen = Diagcolumn[prevcpoint].currentrowindex-
                 Diagcolumn[cpoint].currentrowindex-1;

      evaluateDBcrosspoints(EDtabcolumn,Rtabcolumn,Diagcolumn+cpoint+1,Linear_I,
                          Diagcolumn[cpoint].currentrowindex,coloffset+cpoint+1,
                          useq, Diagcolumn[cpoint].currentrowindex, new_ulen,
                          original_ulen, vseq, vstart + cpoint+1,
                          prevcpoint-cpoint-1, original_vlen,
                          new_left, new_right,
                          matchcost, mismatchcost, gapcost);
    }
    else if (Diagcolumn[prevcpoint].edge == Linear_I)
    {
      dtemp = Diagcolumn[cpoint];
      new_left = MAX(left_dist,
                               -((GtWord)Diagcolumn[prevcpoint].currentrowindex-
                                 (GtWord)Diagcolumn[cpoint].currentrowindex-1));
      new_right = 0;
      new_ulen = Diagcolumn[prevcpoint].currentrowindex-
                 Diagcolumn[cpoint].currentrowindex-1;

      evaluateDBcrosspoints(EDtabcolumn,Rtabcolumn, Diagcolumn+cpoint, Linear_D,
                          Diagcolumn[cpoint].currentrowindex+1,coloffset+cpoint,
                          useq, Diagcolumn[cpoint].currentrowindex+1, new_ulen,
                          original_ulen, vseq, vstart + cpoint,
                          prevcpoint-1-cpoint, original_vlen,
                          new_left, new_right,
                          matchcost, mismatchcost, gapcost);
      Diagcolumn[cpoint]=dtemp;
    }
    else
    {
      /* if (Diagcolumn[cpoint].edge == Linear_X), never reach this line */
      gt_assert(false);
    }
    Diagcolumn[cpoint].lastcpoint  = ctemp;
  }

  /* exception, if first crosspoint != 0-entry */
  if (vstart-coloffset != cpoint)
  {
    if (Diagcolumn[cpoint].edge == Linear_D)
    {
      new_left =  MAX(left_dist,
                -((GtWord)Diagcolumn[cpoint].currentrowindex-(GtWord)ustart-1));
      new_right = MIN(right_dist, (GtWord)cpoint);
      new_ulen = Diagcolumn[cpoint].currentrowindex-ustart-1;

      evaluateDBcrosspoints(EDtabcolumn, Rtabcolumn, Diagcolumn,
                            edge, rowoffset, coloffset,
                            useq, ustart, new_ulen, original_ulen,
                            vseq, vstart, cpoint, original_vlen,
                            new_left, new_right,
                            matchcost, mismatchcost, gapcost);

    }
    else if (Diagcolumn[cpoint].edge == Linear_I)
    {
      new_left = MAX(left_dist,
                 -((GtWord)Diagcolumn[cpoint].currentrowindex-(GtWord)ustart));
      new_right = MIN((GtWord)cpoint-1, right_dist);
      evaluateDBcrosspoints(EDtabcolumn, Rtabcolumn, Diagcolumn,
                            edge, rowoffset, coloffset, useq, ustart,
                            Diagcolumn[cpoint].currentrowindex-ustart,
                            original_ulen, vseq, vstart, cpoint-1,
                            original_vlen, new_left, new_right,
                            matchcost, mismatchcost, gapcost);
    }
    else if (Diagcolumn[cpoint].edge == Linear_R ||
             Diagcolumn[cpoint].edge == Linear_X)
    {
      gt_assert(false); /* this cant be the first crosspoint or
                           it have to be 0-entry */
    }
  }
}

static void gt_calc_diagonalbandalign(const GtUchar *useq,
                                      GtUword ustart, GtUword ulen,
                                      const GtUchar *vseq,
                                      GtUword vstart, GtUword vlen,
                                      GtWord left_dist,
                                      GtWord right_dist,
                                      GtAlignment *align,
                                      GtUword matchcost,
                                      GtUword mismatchcost,
                                      GtUword gapcost)
{
  Diagentry *Diagcolumn;
  GtUword idx, *EDtabcolumn, *Rtabcolumn;

  gt_assert(align != NULL);

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    fprintf(stderr,"ERROR: invalid diagonalband for global alignment "
                   "(ulen: "GT_WU", vlen: "GT_WU")\n"
                   "left_dist <= MIN(0, vlen-ulen) and "
                   "right_dist <= MAX(0, vlen-ulen)\n", ulen, vlen);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  if (ulen == 0UL)
  {
    (void) construct_trivial_insertion_alignment(align,vlen,gapcost);
    return;
  }
  else if (vlen == 0UL)
  {
    (void) construct_trivial_deletion_alignment(align,vlen,gapcost);
    return;
  }
  else if (ulen == 1UL || vlen == 1UL )
  {
    (void) diagonalbandalignment_in_square_space(align, useq, ustart, ulen,
                                                 vseq, vstart, vlen, left_dist,
                                                 right_dist, matchcost,
                                                 mismatchcost, gapcost);
  }
  else
  {
    Diagcolumn = gt_malloc(sizeof *Diagcolumn * (vlen+1));
    EDtabcolumn = gt_malloc(sizeof *EDtabcolumn * (ulen+1));
    Rtabcolumn = gt_malloc(sizeof *Rtabcolumn * (ulen+1));

    /* initialize Diagcolumn */
    for (idx = 0; idx <= vlen; idx++)
    {
      Diagcolumn[idx].lastcpoint = GT_UWORD_MAX;
      Diagcolumn[idx].currentrowindex = GT_UWORD_MAX;
    }

    evaluateDBcrosspoints(EDtabcolumn, Rtabcolumn, Diagcolumn, Linear_X, 0, 0,
                          useq, ustart, ulen, ulen, vseq, vstart, vlen, vlen,
                          left_dist, right_dist,
                          matchcost, mismatchcost, gapcost);

     reconstructalignment_from_Dtab(align,Diagcolumn,ulen, vlen);

    gt_free(Diagcolumn);
    gt_free(EDtabcolumn);
    gt_free(Rtabcolumn);
  }
}

/* compute alignment within a diagonal band */
void gt_computediagonalbandalign(GtAlignment *align,
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
  gt_assert(useq  && vseq);

  /* set new bounds, if left_dist or right_dist is out of sequence */
  left_dist = MAX(-(GtWord) ulen, left_dist);
  right_dist = MIN((GtWord) vlen, right_dist);

  gt_alignment_set_seqs(align,useq+ustart, ulen, vseq+vstart, vlen);
  gt_calc_diagonalbandalign(useq, ustart, ulen, vseq, vstart, vlen,
                            left_dist, right_dist, align,
                            matchcost, mismatchcost, gapcost);
}

void gt_checkdiagonalbandalign(GT_UNUSED bool forward,
                                const GtUchar *useq,
                                GtUword ulen,
                                const GtUchar *vseq,
                                GtUword vlen)
{
  GtUword edist1, edist2, edist3, matchcost = 0, mismatchcost = 1, gapcost = 1;
  GtWord left_dist, right_dist;
  GtAlignment *align;

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

  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);
  gt_calc_diagonalbandalign(useq, 0, ulen, vseq, 0, vlen, left_dist, right_dist,
                            align, matchcost, mismatchcost, gapcost);

  edist3 = gt_alignment_eval_generic_with_score(false, align, matchcost,
                                                mismatchcost, gapcost);

  if (edist2 != edist3)
  {
    fprintf(stderr,"diagonalband_squarespace_distance_only = "GT_WU" != "GT_WU
              " = gt_alignment_eval_with_score\n", edist2, edist3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_alignment_delete(align);
}
