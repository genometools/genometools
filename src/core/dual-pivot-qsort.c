#include <core/types_api.h>
#include <core/assert_api.h>

static void gt_dpqs_exchange(GtUword *input, GtUword idx1,GtUword idx2)
{
  GtUword tmp = input[idx1];

  input[idx1] = input[idx2];
  input[idx2] = tmp;
}

void gt_dual_pivot_quicksort(GtUword *input,GtUword lowindex,GtUword highindex)
{
  GtUword i, lt, gt, pivot1, pivot2;

  gt_assert(lowindex < highindex);
  pivot1 = input[lowindex];
  pivot2 = input[highindex];
  if (pivot1 > pivot2)
  {
    gt_dpqs_exchange(input, lowindex, highindex);
    pivot1 = input[lowindex];
    pivot2 = input[highindex];
    //sort(input, lowindex, highindex);
  } else
  {
    if (pivot1 == pivot2)
    {
      while (pivot1 == pivot2 && lowindex < highindex)
      {
        lowindex++;
        pivot1 = input[lowindex];
      }
    }
  }
  i = lowindex + 1;
  lt = lowindex + 1;
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
  gt_assert(lt > 0);
  gt_dual_pivot_quicksort(input, lowindex, lt - 1);
  gt_assert(gt > 0);
  gt_dual_pivot_quicksort(input, lt + 1, gt - 1);
  gt_dual_pivot_quicksort(input, gt + 1, highindex);
}
