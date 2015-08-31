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
#include <string.h>
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/divmodmul.h"
#include "core/ma_api.h"
#include "extended/affinealign.h"
#include "extended/linearalign_affinegapcost.h"

#include "extended/diagonalbandalign_affinegapcost.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)
typedef struct {
  GtUword Rvalue, Dvalue, Ivalue;
} Atabentry;

typedef struct {
  GtUword lastcpoint, currentrowindex;
  AffineAlignEdge edge;
} Diagnode;

typedef struct {
  Diagnode val_R, val_D, val_I;
} Diagentry;

static GtWord add_safe_max(const GtUword val1, const GtUword val2)
{
  if (val1 != GT_UWORD_MAX && val2 != GT_UWORD_MAX)
  {
    if (val1 > 0 && val2 > 0)
      gt_assert(val1+val2 >= val1 && val1+val2 >= val2);/*check overflow*/
    return val1+val2;
  }
  return GT_UWORD_MAX;
}

/* calculate only distance with diagonalband in square space O(n²) with
 * affine gapcosts */
static GtUword diagonalband_squarespace_affine(const GtUchar *useq,
                                               const GtUword ustart,
                                               const GtUword ulen,
                                               const GtUchar *vseq,
                                               const GtUword vstart,
                                               const GtUword vlen,
                                               const GtWord left_dist,
                                               const GtWord right_dist,
                                               const GtWord matchcost,
                                               const GtWord mismatchcost,
                                               const GtWord gap_opening,
                                               const GtWord gap_extension)
{
  GtUword i,j, low_row, high_row, rcost, Rdist,
          Ddist, Idist, minvalue, distance = GT_UWORD_MAX;
  Atabentry **Atabcolumn;

   if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    return GT_UWORD_MAX;
  }

  low_row = 0;
  high_row = -left_dist;

  Atabcolumn = gt_malloc((sizeof *Atabcolumn)*(ulen+1));
  *Atabcolumn = gt_malloc((sizeof **Atabcolumn)*((vlen+1)*(ulen+1)));
  for (j = 1; j <= ulen; j++)
  {
    Atabcolumn[j] = Atabcolumn[j-1] + vlen + 1;
  }

  /* first column */
  Atabcolumn[0][0].Rvalue = 0;
  Atabcolumn[0][0].Dvalue = gap_opening;
  Atabcolumn[0][0].Ivalue = gap_opening;

  for (i = 1; i <= high_row; i++)
  {
    Atabcolumn[i][0].Rvalue = GT_UWORD_MAX;
    Atabcolumn[i][0].Dvalue = add_safe_max(Atabcolumn[i-1][0].Dvalue,
                                                      gap_extension);
    Atabcolumn[i][0].Ivalue = GT_UWORD_MAX;
  }
  for (; i <= ulen; i++)
  {
    /* invalid values */
    Atabcolumn[i][0].Rvalue = GT_UWORD_MAX;
    Atabcolumn[i][0].Dvalue = GT_UWORD_MAX;
    Atabcolumn[i][0].Ivalue = GT_UWORD_MAX;
  }

  /* next columns */
  for (j = 1; j <= vlen; j++)
  {
    /* below diagonal band*/
    for (i = 0; i <= low_row; i++)
    {
      if (j <= right_dist)
      {
        Rdist = add_safe_max(Atabcolumn[i][j-1].Rvalue,
                             gap_extension + gap_opening);
        Ddist = add_safe_max(Atabcolumn[i][j-1].Dvalue,
                             gap_extension + gap_opening);
        Idist = add_safe_max(Atabcolumn[i][j-1].Ivalue,gap_extension);

        minvalue = MIN3(Rdist, Ddist, Idist);
        Atabcolumn[i][j].Ivalue = minvalue;
        Atabcolumn[i][j].Rvalue = GT_UWORD_MAX;
        Atabcolumn[i][j].Dvalue = GT_UWORD_MAX;
      }
      else{
        Atabcolumn[i][j].Rvalue = GT_UWORD_MAX;
        Atabcolumn[i][j].Dvalue = GT_UWORD_MAX;
        Atabcolumn[i][j].Ivalue = GT_UWORD_MAX;
      }
    }
    if ( j > right_dist)
      low_row++;
    if (high_row < ulen)
      high_row ++;

    /* diagonaldband */
    for (; i <= high_row; i++)
    {
      Rdist = add_safe_max(Atabcolumn[i][j-1].Rvalue,gap_extension+gap_opening);
      Ddist = add_safe_max(Atabcolumn[i][j-1].Dvalue,gap_extension+gap_opening);
      Idist = add_safe_max(Atabcolumn[i][j-1].Ivalue,gap_extension);
      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[i][j].Ivalue = minvalue;

      rcost = useq[ustart+i-1]==vseq[vstart+j-1]? matchcost:mismatchcost;
      Rdist = add_safe_max(Atabcolumn[i-1][j-1].Rvalue, rcost);
      Ddist = add_safe_max(Atabcolumn[i-1][j-1].Dvalue, rcost);
      Idist = add_safe_max(Atabcolumn[i-1][j-1].Ivalue, rcost);
      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[i][j].Rvalue = minvalue;

      Rdist = add_safe_max(Atabcolumn[i-1][j].Rvalue,
                         gap_extension+gap_opening);
      Ddist = add_safe_max(Atabcolumn[i-1][j].Dvalue,gap_extension);
      Idist = add_safe_max(Atabcolumn[i-1][j].Ivalue,
                          gap_extension+gap_opening);
      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[i][j].Dvalue = minvalue;
    }
    /* above diagonal band */
    for (; i <= ulen; i++)
    {
      Atabcolumn[i][j].Rvalue = GT_UWORD_MAX;
      Atabcolumn[i][j].Dvalue = GT_UWORD_MAX;
      Atabcolumn[i][j].Ivalue = GT_UWORD_MAX;
    }
  }
  distance = MIN3(Atabcolumn[ulen][vlen].Rvalue,
                  Atabcolumn[ulen][vlen].Dvalue,
                  Atabcolumn[ulen][vlen].Ivalue);
  gt_free(Atabcolumn[0]);
  gt_free(Atabcolumn);
  return distance;
}

/* calculate only distance with diagonalband in linear space O(n)
 * with affine gapcosts */
static GtUword diagonalband_linear_affine(const GtUchar *useq,
                                          const GtUword ustart,
                                          const GtUword ulen,
                                          const GtUchar *vseq,
                                          const GtUword vstart,
                                          const GtUword vlen,
                                          const GtWord left_dist,
                                          const GtWord right_dist,
                                          const GtWord matchcost,
                                          const GtWord mismatchcost,
                                          const GtWord gap_opening,
                                          const GtWord gap_extension)
{
  GtUword distance, colindex, rowindex, low_row, high_row, width,
  rcost, Rdist,
          Ddist, Idist, minvalue;
  Atabentry *Atabcolumn, Anw, Awe;
  bool last_row = false;

  distance = GT_UWORD_MAX;

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    return GT_UWORD_MAX;
  }

  width = right_dist - left_dist + 1;
  Atabcolumn = gt_malloc(sizeof(*Atabcolumn) * width);

  low_row = 0;
  high_row = -left_dist;
  Atabcolumn[low_row].Rvalue = 0;
  Atabcolumn[low_row].Dvalue = gap_opening;
  Atabcolumn[low_row].Ivalue = gap_opening;

  for (rowindex = low_row+1; rowindex <= high_row; rowindex ++)
  {
    Atabcolumn[rowindex-low_row].Rvalue = GT_UWORD_MAX;
    Atabcolumn[rowindex-low_row].Dvalue = add_safe_max(
                                          Atabcolumn[rowindex-low_row-1].Dvalue,
                                          gap_extension);
    Atabcolumn[rowindex-low_row].Ivalue = GT_UWORD_MAX;
  }
  for (colindex = 1; colindex <= vlen; colindex++)
  {
    Anw = Atabcolumn[0];

    if (colindex > right_dist)
    {
      Awe = Atabcolumn[1];
      low_row++;
    }
    else
      Awe = Atabcolumn[0];
    if (high_row < ulen)
      high_row ++;

    Rdist = add_safe_max(Awe.Rvalue,gap_extension+gap_opening);
    Ddist = add_safe_max(Awe.Dvalue,gap_extension+gap_opening);
    Idist = add_safe_max(Awe.Ivalue,gap_extension);
    minvalue = MIN3(Rdist, Ddist, Idist);
    Atabcolumn[0].Ivalue = minvalue;
    Atabcolumn[0].Rvalue = GT_UWORD_MAX;
    Atabcolumn[0].Dvalue = GT_UWORD_MAX;

    if (low_row > 0 )
    {
      rcost = useq[ustart+rowindex-1] == vseq[vstart+colindex-1]?
              matchcost:mismatchcost;
      Rdist = add_safe_max(Anw.Rvalue, rcost);
      Ddist = add_safe_max(Anw.Dvalue, rcost);
      Idist = add_safe_max(Anw.Ivalue, rcost);
      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[0].Rvalue = minvalue;
    }
    for (rowindex = low_row + 1; rowindex <= high_row; rowindex++)
    {
      Anw = Awe;
      if (!last_row && rowindex == high_row)
      {
        Awe.Rvalue = GT_UWORD_MAX;
        Awe.Dvalue = GT_UWORD_MAX;
        Awe.Ivalue = GT_UWORD_MAX;
      }
      else if (low_row > 0)
        Awe = Atabcolumn[rowindex-low_row+1];
      else
        Awe = Atabcolumn[rowindex-low_row];

      if (rowindex == ulen)
        last_row = true;

      Rdist = add_safe_max(Awe.Rvalue,gap_extension+gap_opening);
      Ddist = add_safe_max(Awe.Dvalue,gap_extension+gap_opening);
      Idist = add_safe_max(Awe.Ivalue,gap_extension);

      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[rowindex-low_row].Ivalue = minvalue;

      rcost = useq[ustart+rowindex-1] == vseq[vstart+colindex-1]?
              matchcost:mismatchcost;
      Rdist = add_safe_max(Anw.Rvalue, rcost);
      Ddist = add_safe_max(Anw.Dvalue, rcost);
      Idist = add_safe_max(Anw.Ivalue, rcost);

      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[rowindex-low_row].Rvalue = minvalue;

      Rdist = add_safe_max(Atabcolumn[rowindex-low_row-1].Rvalue,
                           gap_extension+gap_opening);
      Ddist = add_safe_max(Atabcolumn[rowindex-low_row-1].Dvalue,gap_extension);
      Idist = add_safe_max(Atabcolumn[rowindex-low_row-1].Ivalue,
                           gap_extension+gap_opening);

      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[rowindex-low_row].Dvalue = minvalue;
    }
  }

  distance = MIN3(Atabcolumn[high_row-low_row].Rvalue,
                  Atabcolumn[high_row-low_row].Dvalue,
                  Atabcolumn[high_row-low_row].Ivalue);
  gt_free(Atabcolumn);

  return distance;
}

static void firstEDtabRtabcolumn(Atabentry *Atabcolumn,
                                 Rtabentry *Rtabcolumn,
                                 Diagentry *Diagcolumn,
                                 AffineAlignEdge edge,
                                 const GtWord offset,
                                 GtWord left_dist,
                                 GtWord right_dist,
                                 const GtWord gap_opening,
                                 const GtWord gap_extension)
{
  GtUword rowindex, low_row, high_row;
  GtWord diag;

  diag = GT_DIV2(left_dist + right_dist);
  low_row = 0;
  high_row = -left_dist;
  
  Atabcolumn[low_row].Rvalue = GT_UWORD_MAX;
  Atabcolumn[low_row].Dvalue = GT_UWORD_MAX;
  Atabcolumn[low_row].Ivalue = GT_UWORD_MAX;

  Diagcolumn[0].val_R.currentrowindex = GT_UWORD_MAX;
  Diagcolumn[0].val_R.edge = Affine_X;
  Diagcolumn[0].val_R.lastcpoint = GT_UWORD_MAX;
  
  Diagcolumn[0].val_D.currentrowindex = GT_UWORD_MAX;
  Diagcolumn[0].val_D.edge = Affine_X;
  Diagcolumn[0].val_D.lastcpoint = GT_UWORD_MAX;
  
  Diagcolumn[0].val_I.currentrowindex = GT_UWORD_MAX;
  Diagcolumn[0].val_I.edge = Affine_X;
  Diagcolumn[0].val_I.lastcpoint = GT_UWORD_MAX;
  
  Rtabcolumn[0].val_R.idx = GT_UWORD_MAX;
  Rtabcolumn[0].val_R.edge = Affine_X;
  Rtabcolumn[0].val_D.idx = GT_UWORD_MAX;
  Rtabcolumn[0].val_D.edge = Affine_X;
  Rtabcolumn[0].val_I.idx = GT_UWORD_MAX;
  Rtabcolumn[0].val_I.edge = Affine_X;

  switch (edge) {
  case Affine_R:
    Atabcolumn[low_row].Rvalue = 0;
    if (diag == 0)
    {
      Diagcolumn[0].val_R.currentrowindex = 0 + offset;
      Rtabcolumn[0].val_R.idx = 0;
      Rtabcolumn[0].val_R.edge = Affine_R;
    }
    break;
  case Affine_D:
    Atabcolumn[low_row].Dvalue = 0;
    if (diag == 0)
    {
      Diagcolumn[0].val_D.currentrowindex = 0 + offset;
      Rtabcolumn[0].val_D.idx = 0;
      Rtabcolumn[0].val_D.edge = Affine_D;
    }
    break;
  case Affine_I:
    Atabcolumn[low_row].Ivalue = 0;
    if (diag == 0)
    {
      Diagcolumn[0].val_I.currentrowindex = 0 + offset;
      Rtabcolumn[0].val_I.idx = 0;
      Rtabcolumn[0].val_I.edge = Affine_I;
    }
    break;
  default:
    Atabcolumn[low_row].Rvalue = 0;
    Atabcolumn[low_row].Dvalue = gap_opening;
    Atabcolumn[low_row].Ivalue = gap_opening;
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
    Atabcolumn[rowindex].Rvalue = GT_WORD_MAX;
    Atabcolumn[rowindex].Dvalue = add_safe_max(Atabcolumn[rowindex-1].Dvalue,
                                           gap_extension);
    Atabcolumn[rowindex].Ivalue = GT_WORD_MAX;

    if (diag == -(GtWord)rowindex)
    {      
      Diagcolumn[0].val_D.edge = Affine_X;
      Diagcolumn[0].val_D.lastcpoint = GT_UWORD_MAX;
      Diagcolumn[0].val_D.currentrowindex = rowindex + offset;
      Rtabcolumn[rowindex-low_row].val_D.idx = 0;
      Rtabcolumn[rowindex-low_row].val_D.edge = Affine_D;
    }
    else
    {
      Rtabcolumn[rowindex-low_row] = Rtabcolumn[rowindex-low_row-1];
    }
  }
}

static Rnode evaluateallcolumns(Atabentry *Atabcolumn,
                                Rtabentry *Rtabcolumn,
                                Diagentry *Diagcolumn,
                                AffineAlignEdge edge,
                                const GtWord offset,
                                const GtUchar *useq,
                                const GtUword ustart,
                                const GtUword ulen,
                                const GtUchar *vseq,
                                const GtUword vstart,
                                const GtUword vlen,
                                GtWord left_dist,
                                GtWord right_dist,
                                const GtWord matchcost,
                                const GtWord mismatchcost,
                                const GtWord gap_opening,
                                const GtWord gap_extension)
{
  GtUword colindex, rowindex, Rdist, Ddist, Idist, minvalue, rcost,
          low_row, high_row; /*lowest and highest row between a diagonal band*/
  GtWord diag;
  bool last_row = false;
  Atabentry Anw, Awe;
  Rtabentry Rnw, Rwe;
  Rnode lastcpoint = {GT_UWORD_MAX, Affine_X};

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    gt_assert(false);
  }
  diag = GT_DIV2(left_dist + right_dist);
  printf("diag: "GT_WD"\n",diag);
  low_row = 0;
  high_row = -left_dist;

 /* first column */
  firstEDtabRtabcolumn(Atabcolumn, Rtabcolumn, Diagcolumn, edge, offset,
                     left_dist, right_dist, gap_opening, gap_extension);
  /* next columns */
  for (colindex = 1; colindex <= vlen; colindex++)
  {
    Anw = Atabcolumn[0];
    Rnw = Rtabcolumn[0];

    if (colindex > right_dist)
    {
      Awe = Atabcolumn[1];
      Rwe = Rtabcolumn[1];
      low_row++;
    }
    else
    {
      Awe = Atabcolumn[0];
      Rwe = Rtabcolumn[0];
    }
    if (high_row < ulen)
      high_row ++;

    Rdist = add_safe_max(Atabcolumn[0].Rvalue,gap_extension+gap_opening);
    Ddist = add_safe_max(Atabcolumn[0].Dvalue,gap_extension+gap_opening);
    Idist = add_safe_max(Atabcolumn[0].Ivalue,gap_extension);
  
    minvalue = MIN3(Rdist, Ddist, Idist);
    Atabcolumn[0].Ivalue = minvalue;
    Atabcolumn[0].Rvalue = GT_WORD_MAX;
    Atabcolumn[0].Dvalue = GT_WORD_MAX;
    
    if (diag == (GtWord)colindex - (GtWord)low_row)
    {
      Diagcolumn[colindex].val_R.currentrowindex = GT_UWORD_MAX;
      Diagcolumn[colindex].val_R.edge = Affine_X;
      Diagcolumn[colindex].val_R.lastcpoint = GT_UWORD_MAX;
      
      Diagcolumn[colindex].val_D.currentrowindex = GT_UWORD_MAX;
      Diagcolumn[colindex].val_D.edge = Affine_X;
      Diagcolumn[colindex].val_D.lastcpoint = GT_UWORD_MAX;
      
      Diagcolumn[colindex].val_I.currentrowindex = low_row + offset;
      
      if (minvalue == Rdist)
      {
        Diagcolumn[colindex].val_I.edge = Rwe.val_R.edge;
        Diagcolumn[colindex].val_I.lastcpoint = Rwe.val_R.idx;
      }
      else if (minvalue == Idist)
      {
        Diagcolumn[colindex].val_I.edge = Rwe.val_I.edge;
        Diagcolumn[colindex].val_I.lastcpoint = Rwe.val_I.idx;
      }
      else if (minvalue == Ddist)
      {
        Diagcolumn[colindex].val_I.edge = Rwe.val_D.edge;
        Diagcolumn[colindex].val_I.lastcpoint = Rwe.val_D.idx;
      }
      
      Rtabcolumn[0].val_R.idx = GT_UWORD_MAX;
      Rtabcolumn[0].val_R.edge = Affine_X;
      Rtabcolumn[0].val_D.idx = GT_UWORD_MAX;
      Rtabcolumn[0].val_D.edge = Affine_X;
      Rtabcolumn[0].val_I.idx = colindex;
      Rtabcolumn[0].val_I.edge = Affine_I;
      
    }
    else
    {
      if (minvalue == Rdist)
      {
        Rtabcolumn[0].val_I = Rwe.val_R;
      }
      else if (minvalue == Idist)
      {
       Rtabcolumn[0].val_I = Rwe.val_I;
      }
      else if (minvalue == Ddist)
      {
        Rtabcolumn[0].val_I = Rwe.val_D;
      }
      Rtabcolumn[0].val_D = (Rnode){GT_UWORD_MAX, Affine_X};
      Rtabcolumn[0].val_R = (Rnode){GT_UWORD_MAX, Affine_X};
    }

    /* replacement possible for 0-entry */
    if (low_row > 0 )
    {
      rcost = useq[ustart+low_row-1] == vseq[vstart+colindex-1]? matchcost:mismatchcost;
      Rdist = add_safe_max(Anw.Rvalue, rcost);
      Ddist = add_safe_max(Anw.Dvalue, rcost);
      Idist = add_safe_max(Anw.Ivalue, rcost);

      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[0].Rvalue = minvalue;

      if (diag == (GtWord)colindex - (GtWord)low_row)
      {
        Diagcolumn[colindex].val_R.currentrowindex = low_row + offset;
        
        if (minvalue == Rdist)
        {
          Diagcolumn[colindex].val_R.edge = Rnw.val_R.edge;
          Diagcolumn[colindex].val_R.lastcpoint = Rnw.val_R.idx;
        }
        else if (minvalue == Idist)
        {
          Diagcolumn[colindex].val_R.edge = Rnw.val_I.edge;
          Diagcolumn[colindex].val_R.lastcpoint = Rnw.val_I.idx;
        }
        else if (minvalue == Ddist)
        {
          Diagcolumn[colindex].val_R.edge = Rnw.val_D.edge;
          Diagcolumn[colindex].val_R.lastcpoint = Rnw.val_D.idx;
        }

        Rtabcolumn[0].val_R.idx = colindex;
        Rtabcolumn[0].val_R.edge = Affine_R;
      }
      else
      {
        if (minvalue == Rdist)
          Rtabcolumn[0].val_R = Rnw.val_R;
        else if (minvalue == Idist)
           Rtabcolumn[0].val_R = Rnw.val_I;
        else if (minvalue == Ddist)
           Rtabcolumn[0].val_R = Rnw.val_D;
      }
    }
    for (rowindex = low_row + 1; rowindex <= high_row; rowindex++)
    {
      Anw = Awe;
      Rnw = Rwe;

      if (!last_row && rowindex == high_row)
      {/* prev is outside of diagonalband*/
        Awe = (Atabentry){GT_UWORD_MAX,GT_UWORD_MAX,GT_UWORD_MAX};
        Rwe.val_R = (Rnode){GT_UWORD_MAX,Affine_X};
        Rwe.val_D = (Rnode){GT_UWORD_MAX,Affine_X};
        Rwe.val_I = (Rnode){GT_UWORD_MAX,Affine_X};
      }
      else if (low_row > 0)
      {/* shifted diagonalband*/
        Awe = Atabcolumn[rowindex-low_row+1];
        Rwe = Rtabcolumn[rowindex-low_row+1];
      }
      else
      {/* normaly prev*/
        Awe = Atabcolumn[rowindex-low_row];
        Rwe = Rtabcolumn[rowindex-low_row];
      }
      if (rowindex == ulen)
        last_row = true;
      /* insertion */
      Rdist = add_safe_max(Awe.Rvalue,gap_extension+gap_opening);
      Ddist = add_safe_max(Awe.Dvalue,gap_extension+gap_opening);
      Idist = add_safe_max(Awe.Ivalue,gap_extension);

      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[rowindex-low_row].Ivalue = minvalue;
      if (diag == (GtWord)colindex - (GtWord)rowindex)
      {
        if (minvalue == Rdist)
        {
          Diagcolumn[colindex].val_I.edge = Rwe.val_R.edge;
          Diagcolumn[colindex].val_I.lastcpoint = Rwe.val_R.idx;
        }
        else if (minvalue == Idist)
        {
          Diagcolumn[colindex].val_I.edge = Rwe.val_I.edge;
          Diagcolumn[colindex].val_I.lastcpoint = Rwe.val_I.idx;
        }
        else if (minvalue == Ddist)
        {
          Diagcolumn[colindex].val_I.edge = Rwe.val_D.edge;
          Diagcolumn[colindex].val_I.lastcpoint = Rwe.val_D.idx;
        }
        Diagcolumn[colindex].val_I.currentrowindex = rowindex+offset;
        Rtabcolumn[rowindex-low_row].val_I.idx = colindex;
        Rtabcolumn[rowindex-low_row].val_I.edge = Affine_I;
      }
      else
      {
        if (minvalue == Rdist)
          Rtabcolumn[rowindex-low_row].val_I = Rwe.val_R;
        else if (minvalue == Idist)
           Rtabcolumn[rowindex-low_row].val_I = Rwe.val_I;
        else if (minvalue == Ddist)
           Rtabcolumn[rowindex-low_row].val_I = Rwe.val_D;
      }
      /* replacement */
      rcost = useq[ustart+rowindex-1]==vseq[vstart+colindex-1]? matchcost:mismatchcost;
      Rdist = add_safe_max(Anw.Rvalue, rcost);
      Ddist = add_safe_max(Anw.Dvalue, rcost);
      Idist = add_safe_max(Anw.Ivalue, rcost);
      minvalue = MIN3(Rdist, Ddist, Idist);
    
      Atabcolumn[rowindex-low_row].Rvalue = minvalue;
      if (diag == (GtWord)colindex - (GtWord)rowindex)
      {
        if (minvalue == Rdist)
        {
          Diagcolumn[colindex].val_R.edge = Rnw.val_R.edge;
          Diagcolumn[colindex].val_R.lastcpoint = Rnw.val_R.idx;
        }
        else if (minvalue == Idist)
        {
          Diagcolumn[colindex].val_R.edge = Rnw.val_I.edge;
          Diagcolumn[colindex].val_R.lastcpoint = Rnw.val_I.idx;
        }
        else if (minvalue == Ddist)
        {
          Diagcolumn[colindex].val_R.edge = Rnw.val_D.edge;
          Diagcolumn[colindex].val_R.lastcpoint = Rnw.val_D.idx;
        }

        Diagcolumn[colindex].val_R.currentrowindex = rowindex+offset;
        Rtabcolumn[rowindex-low_row].val_R.idx = colindex;
        Rtabcolumn[rowindex-low_row].val_R.edge = Affine_R;
      }
      else
      {
        if (minvalue == Rdist)
          Rtabcolumn[rowindex-low_row].val_R = Rnw.val_R;
        else if (minvalue == Idist)
           Rtabcolumn[rowindex-low_row].val_R = Rnw.val_I;
        else if (minvalue == Ddist)
           Rtabcolumn[rowindex-low_row].val_R = Rnw.val_D;
      }
      
      /* deletion */
      Rdist = add_safe_max(Atabcolumn[rowindex-low_row-1].Rvalue,
                           gap_extension+gap_opening);
      Ddist = add_safe_max(Atabcolumn[rowindex-low_row-1].Dvalue,gap_extension);
      Idist = add_safe_max(Atabcolumn[rowindex-low_row-1].Ivalue,
                           gap_extension+gap_opening);
    
      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[rowindex].Dvalue = minvalue;

      if (diag == (GtWord)colindex - (GtWord)rowindex)
      {
        if (minvalue == Rdist)
        {
          Diagcolumn[colindex].val_D.edge = Rtabcolumn[rowindex-low_row-1].val_R.edge;
          Diagcolumn[colindex].val_D.lastcpoint = Rtabcolumn[rowindex-low_row-1].val_R.idx;
        }
        else if (minvalue == Idist)
        {
          Diagcolumn[colindex].val_D.edge = Rtabcolumn[rowindex-low_row-1].val_I.edge;
          Diagcolumn[colindex].val_D.lastcpoint = Rtabcolumn[rowindex-low_row-1].val_I.idx;
        }
        else if (minvalue == Ddist)
        {
          Diagcolumn[colindex].val_D.edge = Rtabcolumn[rowindex-low_row-1].val_D.edge;
          Diagcolumn[colindex].val_D.lastcpoint = Rtabcolumn[rowindex-low_row-1].val_D.idx;
        }

        Diagcolumn[colindex].val_D.currentrowindex = rowindex+offset;
        Rtabcolumn[rowindex-low_row].val_D.idx = colindex;
        Rtabcolumn[rowindex-low_row].val_D.edge = Affine_D;
      }
      else
      {
        if (minvalue == Rdist)
          Rtabcolumn[rowindex-low_row].val_D = Rtabcolumn[rowindex-low_row-1].val_R;
        else if (minvalue == Idist)
           Rtabcolumn[rowindex-low_row].val_D = Rtabcolumn[rowindex-low_row-1].val_I;
        else if (minvalue == Ddist)
           Rtabcolumn[rowindex-low_row].val_D = Rtabcolumn[rowindex-low_row-1].val_D;
      }
      
    }
  }
  /* last crosspoint of optimal path */
  Rdist = Atabcolumn[high_row-low_row].Rvalue;
  Ddist = Atabcolumn[high_row-low_row].Dvalue;
  Idist = Atabcolumn[high_row-low_row].Ivalue;
  minvalue = MIN3(Rdist, Ddist, Idist);

  if (minvalue == Rdist)
    lastcpoint = Rtabcolumn[high_row-low_row].val_R;
  else if (minvalue == Idist)
    lastcpoint = Rtabcolumn[high_row-low_row].val_I;
  else if (minvalue == Ddist)
    lastcpoint = Rtabcolumn[high_row-low_row].val_D;
  
  return lastcpoint;
}

static void evaluatecrosspoints(Atabentry *Atabcolumn,
                                Rtabentry *Rtabcolumn,
                                Diagentry *Diagcolumn,
                                AffineAlignEdge edge,
                                const GtUword rowoffset,
                                GT_UNUSED const GtUword coloffset,//unsused
                                const GtUchar *useq,
                                const GtUword ustart,
                                const GtUword ulen,
                                const GtUchar *vseq,
                                const GtUword vstart,
                                const GtUword vlen,
                                const GtWord left_dist,
                                const GtWord right_dist,
                                const GtWord matchcost,
                                const GtWord mismatchcost,
                                const GtWord gap_opening,
                                const GtWord gap_extension)
{//TODO recursive
  Rnode cpoint;

  cpoint = evaluateallcolumns(Atabcolumn, Rtabcolumn, Diagcolumn, edge,
                              rowoffset, useq, ustart, ulen, vseq, vstart, vlen,
                              left_dist, right_dist,
                              matchcost, mismatchcost, gap_opening, gap_extension);
  printf("idx: "GT_WU", edge: %d\n", cpoint.idx, cpoint.edge);
}

static void gt_calc_diagonalbandaffinealign(const GtUchar *useq,
                                            GtUword ustart, GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart, GtUword vlen,
                                            GtWord left_dist,
                                            GtWord right_dist,
                                            GtAlignment *align,
                                            const GtWord matchcost,
                                            const GtWord mismatchcost,
                                            const GtWord gap_opening,
                                            const GtWord gap_extension)
{
  Atabentry *Atabcolumn;
  Rtabentry *Rtabcolumn;
  Diagentry *Diagcolumn;
  GtUword idx;

  gt_assert(align != NULL);

  if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    fprintf(stderr,"invalid diagonalband for global alignment\n");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  Diagcolumn = gt_malloc(sizeof *Diagcolumn * (vlen+1));
  Atabcolumn = gt_malloc(sizeof *Atabcolumn * (ulen+1));
  Rtabcolumn = gt_malloc(sizeof *Rtabcolumn * (ulen+1));

  /* initialize Diagcolumn */
  for (idx = 0; idx <= vlen; idx++)
  {
    Diagcolumn[idx].val_R = (Diagnode){GT_UWORD_MAX, GT_UWORD_MAX, Affine_X};
    Diagcolumn[idx].val_D = (Diagnode){GT_UWORD_MAX, GT_UWORD_MAX, Affine_X};
    Diagcolumn[idx].val_I = (Diagnode){GT_UWORD_MAX, GT_UWORD_MAX, Affine_X};
  }

  evaluatecrosspoints(Atabcolumn, Rtabcolumn, Diagcolumn, Affine_X, 0, 0,
                      useq, ustart, ulen, vseq, vstart, vlen,
                      left_dist, right_dist,
                      matchcost, mismatchcost, gap_opening,gap_extension);
  /* reconstruct alignment */
  //TODO

  gt_free(Diagcolumn);
  gt_free(Atabcolumn);
  gt_free(Rtabcolumn);
}

void gt_computediagnoalbandaffinealign(GtAlignment * align,
                                       const GtUchar *useq,
                                       GtUword ustart, GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart, GtUword vlen,
                                       GtWord left_dist,
                                       GtWord right_dist,
                                       const GtWord matchcost,
                                       const GtWord mismatchcost,
                                       const GtWord gap_opening,
                                       const GtWord gap_extension)
{

  gt_assert(useq  && vseq);
  if (matchcost < 0 || mismatchcost < 0 || gap_opening < 0 || gap_extension < 0)
  {
    fprintf(stderr,"invalid cost value\n");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  /* set new bounds, if left_dist or right_dist is out of sequence */
  left_dist = MAX(-(GtWord) ulen,left_dist);
  right_dist = MIN((GtWord) vlen,right_dist);

  align = gt_alignment_new_with_seqs(useq + ustart, ulen, vseq + vstart, vlen);
  gt_calc_diagonalbandaffinealign(useq, ustart, ulen, vseq, vstart, vlen,
                            left_dist, right_dist, align,
                            matchcost, mismatchcost, gap_opening, gap_extension);
}

void gt_checkdiagnonalbandaffinealign(GT_UNUSED bool forward,
                                const GtUchar *useq,
                                GtUword ulen,
                                const GtUchar *vseq,
                                GtUword vlen)
{
  GtUword affine_cost1, affine_cost2, affine_cost3;
  GtWord left_dist, right_dist, matchcost = 0, mismatchcost = 4,
         gap_opening = 4, gap_extension = 1;
  GtAlignment *align_square;

  if (strchr((const char*)useq, LINEAR_EDIST_GAP))
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (strchr((const char*)vseq, LINEAR_EDIST_GAP))
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  /* set left and right to set diagonalband to the whole matrix */
  left_dist = -ulen;
  right_dist = vlen;
  affine_cost1 = diagonalband_squarespace_affine(useq, 0, ulen, vseq, 0, vlen,
                                                 left_dist, right_dist,
                                                 matchcost, mismatchcost,
                                                 gap_opening, gap_extension);
  align_square = gt_affinealign(useq, ulen, vseq, vlen, matchcost,
                                mismatchcost, gap_opening, gap_extension);
  affine_cost2 = gt_alignment_eval_with_affine_score(align_square, matchcost,
                                                     mismatchcost, gap_opening,
                                                     gap_extension);
  gt_alignment_delete(align_square);

  if (affine_cost1 != affine_cost2)
  {
    fprintf(stderr,"diagonalband_squarespace_affine = "GT_WU
            " != "GT_WU" = gt_affinealign\n", affine_cost1, affine_cost2);

    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  affine_cost3 = diagonalband_linear_affine(useq, 0, ulen, vseq, 0, vlen,
                                            left_dist, right_dist,
                                            matchcost, mismatchcost,
                                            gap_opening, gap_extension);
  if (affine_cost3 != affine_cost2)
  {
    fprintf(stderr,"diagonalband_linear_affine = "GT_WU
            " != "GT_WU" = gt_affinealign\n", affine_cost3, affine_cost2);

    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
