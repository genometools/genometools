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
#include "extended/diagonalband_align.h"


#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)
typedef enum {
  R,
  D,
  I,
  X /*unknown*/
} Edge;

typedef struct {
  GtWord Rvalue, Dvalue, Ivalue, totalvalue;

  Edge Redge,
       Dedge,
       Iedge;
}Atabentry;

static GtWord add_safe_max(const GtWord val1, const GtWord val2)
{
  if (val1 != GT_WORD_MAX && val2 != GT_WORD_MAX)
  {
     //if (val1 > 0 && val2 > 0)
       //gt_assert(val1+val2 >= val1 && val1+val2 >= val2);/*check overflow*/
     return val1+val2;
  }

    return GT_WORD_MAX;
}

/* affine distance_only */
GtUword diagonalband_affine_distance_only(const GtUchar *useq,
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

  for(rowindex = low_row+1; rowindex <= high_row; rowindex ++)
  {
    Atabcolumn[rowindex-low_row].Rvalue = GT_WORD_MAX;
    Atabcolumn[rowindex-low_row].Dvalue = 
             add_safe_max(Atabcolumn[rowindex-low_row-1].Dvalue, gap_extension);
    Atabcolumn[rowindex-low_row].Ivalue = GT_WORD_MAX;
  }

  for(colindex = 1; colindex <= vlen; colindex++)
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
      rcost = (useq[ustart+low_row-1] == vseq[vstart+colindex-1])? matchcost:mismatchcost;
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

    for(rowindex = low_row + 1; rowindex <= high_row; rowindex++)
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
      if (rowindex != high_row)
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
}

/*void gt_checkaffinelinearspace_local(GT_UNUSED bool forward,
                                       const GtUchar *useq,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vlen)
{

}*/
