/*
  Copyright (c) 2014 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

/*
  The code was adapted from the Java code in
  http://algs4.cs.princeton.edu/23quicksort/QuickDualPivot.java.html
  Can be optimized by determining better pivot elements and
  using insertion sort for small subarrays to soet.
*/

#include "core/types_api.h"
#include "core/assert_api.h"
#include "core/dual-pivot-qsort.h"

static void gt_dpqs_exchange(GtUword *input, GtUword idx1,GtUword idx2)
{
  GtUword tmp = input[idx1];

  input[idx1] = input[idx2];
  input[idx2] = tmp;
}

static void gt_rec_dual_pivot_quicksort(GtUword *input,GtUword lowindex,
                                        GtUword highindex)
{
  GtUword i, lt, gt, pivot1, pivot2;

  if (lowindex >= highindex)
  {
    return;
  }
  if (input[lowindex] > input[highindex])
  {
    gt_dpqs_exchange(input, lowindex, highindex);
  }
  pivot1 = input[lowindex];
  pivot2 = input[highindex];
  gt_assert(pivot1 <= pivot2);
  i = lowindex + 1;
  lt = lowindex + 1;
  gt_assert(highindex > 0);
  gt = highindex - 1;
  while (i <= gt)
  {
    if (input[i] < pivot1)
    {
      gt_dpqs_exchange(input, i++, lt++);
    } else
    {
      if (pivot2 < input[i])
      {
        gt_dpqs_exchange(input, i, gt--);
      } else
      {
        i++;
      }
    }
  }
  gt_assert(lt > 0);
  gt_dpqs_exchange(input, lowindex, --lt);
  gt_dpqs_exchange(input, highindex, ++gt);
  if (lt >= 1UL)
  {
    gt_rec_dual_pivot_quicksort(input, lowindex, lt - 1);
  }
  if (gt >= 1UL && input[lt] < input[gt])
  {
    gt_rec_dual_pivot_quicksort(input, lt + 1, gt - 1);
  }
  gt_rec_dual_pivot_quicksort(input, gt + 1, highindex);
}

void gt_dual_pivot_qsort(GtUword *source,GtUword len)
{
  if (len >= 2UL)
  {
    gt_rec_dual_pivot_quicksort(source,0,len-1);
  }
}
