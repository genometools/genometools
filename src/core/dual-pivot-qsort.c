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
  using insertion sort for small subarrays to sort.
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

static void gt_dual_insertion_sort(GtUword *input,GtUword lowindex,
                                   GtUword highindex)
{
  GtUword pm;

  for (pm = lowindex + 1; pm <= highindex; pm++)
  {
    GtUword pl;

    for (pl = pm; pl > lowindex && input[pl-1] > input[pl]; pl--)
    {
      gt_dpqs_exchange(input, pl, pl-1);
    }
  }
}

static void gt_rec_dual_pivot_quicksort(GtUword *input,GtUword lowindex,
                                        GtUword highindex)
{
  if (lowindex >= highindex)
  {
    return;
  } else
  {
    GtUword len = highindex - lowindex + 1;

    if (len < 5UL)
    {
      gt_dual_insertion_sort(input,lowindex,highindex);
    } else
    {
      GtUword i, lt, gt, pivot1, pivot2;

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
  }
}

static void gt_dual_pivots_get(GtUword *pivot1,GtUword *pivot2,
                               GtUword *input,GtUword lowindex,
                               GtUword highindex)
{
  GtUword s7, idx, len = highindex - lowindex + 1, e[6], values[5];

  s7 = (len >> 3) + (len >> 6) + 1;
  e[3] = (lowindex + highindex) >> 1;
  gt_assert(e[3] >= s7);
  e[2] = e[3] - s7;
  e[1] = e[2] - s7;
  if (len == 8UL)
  {
    e[1]++;
  }
  e[4] = e[3] + s7;
  e[5] = e[4] + s7;
  for (idx = (GtUword) 1; idx <= (GtUword) 5; idx++)
  {
    gt_assert(lowindex <= e[idx] && e[idx] <= highindex);
    values[idx-1] = input[e[idx]];
  }
  gt_dual_insertion_sort(values,0,(GtUword) 4);
  for (idx = (GtUword) 1; idx <= (GtUword) 5; idx++)
  {
    input[e[idx]] = values[idx-1];
  }
  *pivot1 = input[e[2]];
  *pivot2 = input[e[4]];
  input[lowindex] = *pivot1;
  input[highindex] = *pivot2;
}

/* This does not work at the moment. Needs to be fixed */

void gt_rec_dual_pivot_quicksort_fast(GtUword *input,GtUword lowindex,
                                      GtUword highindex)
{
  if (lowindex >= highindex)
  {
    return;
  } else
  {
    GtUword len = highindex - lowindex + 1;

    if (len <= 6UL)
    {
      gt_dual_insertion_sort(input,lowindex,highindex);
    } else
    {
      GtUword idx, lt, gt, pivot1, pivot2;

      gt_dual_pivots_get(&pivot1,&pivot2,input,lowindex,highindex);
      gt_assert(pivot1 <= pivot2);
      lt = lowindex + 1;
      gt_assert(highindex > 0);
      gt = highindex - 1;
      for (idx = lowindex + 1; idx <= gt; idx++)
      {
        GtUword currentvalue = input[idx];

        if (currentvalue < pivot1)
        {
          input[idx] = input[lt];
          input[lt++] = currentvalue;
        } else
        {
          if (currentvalue >= pivot2)
          {
            while (input[gt] > pivot2 && idx < gt)
            {
              gt--;
            }
            if (input[gt] < pivot1)
            {
              input[idx] = input[lt];
              input[lt++] = input[gt];
            } else
            {
              input[idx] = input[gt];
            }
            input[gt] = currentvalue;
            gt--;
          }
        }
      }
      gt_assert(lt > 0);
      lt--;
      gt++;
      input[lowindex] = input[lt];
      input[lt] = pivot1;
      input[highindex] = input[gt];
      input[gt] = pivot2;
      if (lt >= 1UL)
      {
        gt_rec_dual_pivot_quicksort_fast(input, lowindex, lt - 1);
      }
      gt_rec_dual_pivot_quicksort_fast(input, lt, gt);
      gt_rec_dual_pivot_quicksort_fast(input, gt + 1, highindex);
    }
  }
}

void gt_dual_pivot_qsort(GtUword *source,GtUword len)
{
  if (len >= 2UL)
  {
    gt_rec_dual_pivot_quicksort(source,0,len-1);
  }
}
