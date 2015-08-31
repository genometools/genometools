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

#include "extended/diagonalbandalign_affinegapcost.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)
typedef struct {
  GtUword Rvalue, Dvalue, Ivalue;
} Atabentry;

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

void gt_checkdiagnonalbandaffinealign(GT_UNUSED bool forward,
                                const GtUchar *useq,
                                GtUword ulen,
                                const GtUchar *vseq,
                                GtUword vlen)
{
  GtUword affine_cost1, affine_cost2, affine_cost3;
  GtWord left_dist, right_dist;
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
                                                 0, 4, 4, 1);
  align_square = gt_affinealign(useq, ulen, vseq, vlen, 0, 4, 4, 1);
  affine_cost2 = gt_alignment_eval_with_affine_score(align_square, 0, 4, 4, 1);
  gt_alignment_delete(align_square);

  if (affine_cost1 != affine_cost2)
  {
    fprintf(stderr,"diagonalband_squarespace_affine = "GT_WU
            " != "GT_WU" = gt_affinealign\n", affine_cost1, affine_cost2);

    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  affine_cost3 = diagonalband_linear_affine(useq, 0, ulen, vseq, 0, vlen,
                                            left_dist, right_dist,
                                             0, 4, 4, 1);
  if (affine_cost3 != affine_cost2)
  {
    fprintf(stderr,"diagonalband_linear_affine = "GT_WU
            " != "GT_WU" = gt_affinealign\n", affine_cost3, affine_cost2);

    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
