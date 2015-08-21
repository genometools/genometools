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
#include "core/assert_api.h"
#include "core/minmax.h"
#include "core/error.h"
#include "core/types_api.h"
#include "core/divmodmul.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "match/squarededist.h"
#include "extended/diagonalband_align.h"
#include "extended/linearalign.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

static GtUword add_safe_max(const GtUword val1, const GtUword val2)
{
  if (val1 != GT_UWORD_MAX && val2 != GT_UWORD_MAX)
  {
     if (val1 > 0 && val2 > 0)
       gt_assert(val1+val2 >= val1 && val1+val2 >= val2);/*check overflow*/
     return val1+val2;
  }

    return GT_UWORD_MAX;
}

/* calculate only distance with diagonalband in sqaure space O(n²) */
static GtUword diagonalband_squarespace_distance_only(const GtUchar *useq,
                                                      const GtUword ustart,
                                                      const GtUword ulen,
                                                      const GtUchar *vseq,
                                                      const GtUword vstart,
                                                      const GtUword vlen,
                                                      const GtWord left_dist,
                                                      const GtWord right_dist,
                                                      const GtWord matchcost,
                                                      const GtWord mismatchcost,
                                                      const GtWord gapcost)
{
  GtUword **E, i,j, val, low_row, high_row, distance = GT_UWORD_MAX;

   if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
      (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
  {
    return GT_UWORD_MAX;
  }

  low_row = 0;
  high_row = -left_dist;

  E = gt_malloc((sizeof **E)*(ulen+1));
  *E = gt_malloc((sizeof *E)*((vlen+1)*(ulen+1)));
  for (j = 1; j <= ulen; j++)
  {
    E[j] = E[j-1]+vlen+1;
  }

  E[0][0] = 0;
  for (i = 1; i <= high_row; i++)
  {
      E[i][0] = add_safe_max(E[i-1][0], gapcost);
  }
  for (; i <= ulen; i++)
  {
      E[i][0] = GT_UWORD_MAX;
  }

  for (j = 1; j <= vlen; j++)
  {
    for (i = 0; i <= low_row; i++)
    {
      if (j < right_dist)
      {
        E[i][j] = add_safe_max(E[i][j-1], gapcost);
      }
      else{
        E[i][j] = GT_UWORD_MAX;
      }
    }
    if ( j > right_dist)
      low_row++;
    if (high_row < ulen)
      high_row ++;
    for (; i <= high_row; i++)
    {
      E[i][j] = add_safe_max(E[i][j-1], gapcost);

      if ((val = add_safe_max(E[i-1][j-1],(useq[ustart+i-1] == vseq[vstart+j-1]?
                                matchcost : mismatchcost)))
          <= E[i][j])
      {
        E[i][j] = val;
      }

      if ((val = add_safe_max(E[i-1][j],gapcost)) < E[i][j])
      {
        E[i][j] = val;
      }
    }
    for (; i <= ulen; i++)
      E[i][j] = GT_UWORD_MAX;
  }
  distance = E[ulen][vlen];

  gt_free(E[0]);
  gt_free(E);
  return distance;
}

/* calculate only distance with diagonalband in linear space O(n) */
GtUword diagonalband_linear_distance_only(const GtUchar *useq,
                                          const GtUword ustart,
                                          const GtUword ulen,
                                          const GtUchar *vseq,
                                          const GtUword vstart,
                                          const GtUword vlen,
                                          const GtWord left_dist,
                                          const GtWord right_dist,
                                          const GtWord matchcost,
                                          const GtWord mismatchcost,
                                          const GtWord gapcost)
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
    EDtabcolumn[0] = add_safe_max(westEDtabentry, 1);

    if (low_row > 0 )
    {
      if ((val = add_safe_max(northwestEDtabentry,
                             (useq[ustart+low_row-1] == vseq[vstart+colindex-1]?
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
      EDtabcolumn[rowindex-low_row] = add_safe_max(westEDtabentry, gapcost);

      val = add_safe_max(northwestEDtabentry,
                        (useq[ustart+rowindex-1] == vseq[vstart+colindex-1] ?
                         matchcost : mismatchcost));
      if (val <= EDtabcolumn[rowindex-low_row])
        EDtabcolumn[rowindex-low_row] = val;

      if ((val = add_safe_max(EDtabcolumn[rowindex-low_row-1], gapcost))
                                       <= EDtabcolumn[rowindex-low_row])
        EDtabcolumn[rowindex-low_row] = val;
    }
  }

  distance = EDtabcolumn[high_row-low_row];
  gt_free(EDtabcolumn);

  return distance;
}

/* affine distance_only */
/*GtWord diagonalband_affine_distance_only(const GtUchar *useq,
                                          const GtUword ustart,
                                          const GtUword ulen,
                                          const GtUchar *vseq,
                                          const GtUword vstart,
                                          const  GtUword vlen,
                                          const GtWord left_dist,
                                          const GtWord right_dist,
                                          const GtWord matchcost,
                                          const GtWord mismatchcost,
                                          const GtWord gap_opening,
                                          const GtWord gap_extension)
{
  GtWord distance, colindex, rowindex, low_row, high_row, width,
         rcost, Rdist, Ddist, Idist, minvalue, maxvalue;
  Atabentry *Atabcolumn, northwestAtabentry, westAtabentry;

  minvalue = MIN(0, vlen-ulen);
  gt_assert(left_dist <= minvalue);
  maxvalue =  MAX(0, vlen-ulen);
  gt_assert(right_dist >= maxvalue);
  distance = GT_WORD_MAX;
  width = right_dist - left_dist + 1;
  Atabcolumn = gt_malloc(sizeof(*Atabcolumn) * width);

  low_row = 0;
  high_row = -left_dist;
  Atabcolumn[low_row].Rvalue = 0;
  Atabcolumn[low_row].Dvalue = gap_opening;
  Atabcolumn[low_row].Ivalue = gap_opening;
  Atabcolumn[high_row+1].Rvalue = GT_WORD_MAX;
  Atabcolumn[high_row+1].Dvalue = GT_WORD_MAX;
  Atabcolumn[high_row+1].Ivalue = GT_WORD_MAX;

  for (rowindex = low_row+1; rowindex <= high_row; rowindex ++)
  {
    Atabcolumn[rowindex-low_row].Rvalue = GT_WORD_MAX;
    Atabcolumn[rowindex-low_row].Dvalue =
             add_safe_max(Atabcolumn[rowindex-low_row-1].Dvalue, gap_extension);
    Atabcolumn[rowindex-low_row].Ivalue = GT_WORD_MAX;
  }

  for (colindex = 1; colindex <= vlen; colindex++)
  {
    northwestAtabentry = Atabcolumn[0];
    westAtabentry = Atabcolumn[0];
    if (colindex > right_dist)
    {
      westAtabentry = Atabcolumn[1];
      low_row++;
    }
    if (high_row < ulen)
      high_row ++;

    if (low_row > 0)
    {
      rcost = (useq[ustart+low_row-1] == vseq[vstart+colindex-1])?
               matchcost:mismatchcost;
      Rdist = add_safe_max(northwestAtabentry.Rvalue,rcost);
      Ddist = add_safe_max(northwestAtabentry.Dvalue,rcost);
      Idist = add_safe_max(northwestAtabentry.Ivalue,rcost);

      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[0].Rvalue = minvalue;
    }
    else
      Atabcolumn[0].Rvalue = GT_WORD_MAX;

    Rdist = add_safe_max(westAtabentry.Rvalue,gap_extension+gap_opening);
    Ddist = add_safe_max(westAtabentry.Dvalue,gap_extension+gap_opening);
    Idist = add_safe_max(westAtabentry.Ivalue,gap_extension);

    minvalue = MIN3(Rdist, Ddist, Idist);
    Atabcolumn[0].Dvalue = GT_WORD_MAX;
    Atabcolumn[0].Ivalue = minvalue;

    for (rowindex = low_row + 1; rowindex <= high_row; rowindex++)
    {
      northwestAtabentry = westAtabentry;
      if (low_row > 0)
        westAtabentry = Atabcolumn[rowindex-low_row + 1];
      else
       westAtabentry = Atabcolumn[rowindex-low_row];

      rcost = (useq[ustart+rowindex-1] == vseq[vstart+colindex-1])?
              matchcost:mismatchcost;
      Rdist = add_safe_max(northwestAtabentry.Rvalue, rcost);
      Ddist = add_safe_max(northwestAtabentry.Dvalue, rcost);
      Idist = add_safe_max(northwestAtabentry.Ivalue, rcost);

      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[rowindex-low_row].Rvalue = minvalue;

      Rdist = add_safe_max(Atabcolumn[rowindex-low_row-1].Rvalue,
                         gap_extension+gap_opening);
      Ddist = add_safe_max(Atabcolumn[rowindex-low_row-1].Dvalue,gap_extension);
      Idist = add_safe_max(Atabcolumn[rowindex-low_row-1].Ivalue,
                        gap_extension+gap_opening);

      minvalue = MIN3(Rdist, Ddist, Idist);
      Atabcolumn[rowindex-low_row].Dvalue = minvalue;

      Rdist = add_safe_max(westAtabentry.Rvalue,gap_extension+gap_opening);
      Ddist = add_safe_max(westAtabentry.Dvalue,gap_extension+gap_opening);
      if (rowindex != high_row//|| rowindex ==ulen
      {
        Idist = add_safe_max(westAtabentry.Ivalue,gap_extension);
        minvalue = MIN3(Rdist, Ddist, Idist);
        Atabcolumn[rowindex-low_row].Ivalue = minvalue;
      }
      else
      {
        Atabcolumn[rowindex-low_row].Ivalue = GT_WORD_MAX;
      }
    }
  }
  distance = MIN3(Atabcolumn[high_row-low_row].Rvalue,
                  Atabcolumn[high_row-low_row].Dvalue,
                  Atabcolumn[high_row-low_row].Ivalue);
  gt_free(Atabcolumn);
  return distance;
}*/

void gt_checkdiagnonalbandalign(GT_UNUSED bool forward,
                                       const GtUchar *useq,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vlen)
{
  GtUword edist1, edist2;
  GtWord left_dist, right_dist;

  /* example to set left_dist and right_dist */
  if ((GtWord)vlen-(GtWord)ulen > 0)
    left_dist = 0;
  else
    left_dist = (GtWord)vlen-(GtWord)ulen-1;
  if ((GtWord)vlen-(GtWord)ulen > 0)
    right_dist = (GtWord)vlen-(GtWord)ulen+2;
  else
    right_dist = 0;

  edist1 = diagonalband_linear_distance_only(useq, 0, ulen,
                                   vseq,0, vlen, left_dist, right_dist,
                                   0,1,1);
  edist2 = diagonalband_squarespace_distance_only(useq, 0, ulen,
                                   vseq,0, vlen, left_dist, right_dist,
                                   0,1,1);

  if (edist1 != edist2)
  {
    fprintf(stderr,"diagonalband_linear_distance_only = "GT_WU" != "GT_WU
              " = diagonalband_squarespace_distance_only\n", edist1, edist2);

    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

}
