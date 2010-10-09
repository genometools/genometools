/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "core/divmodmul.h"
#include "eis-voiditf.h"
#include "pckdfs.h"

typedef struct
{
  unsigned long lowerbound,
                upperbound;
  unsigned int depth;
  bool endwithspecial;
} Dfs_Boundsatdepth;

GT_DECLAREARRAYSTRUCT(Dfs_Boundsatdepth);

static unsigned long dfsnonspecialwidth(const GtArrayBoundswithchar *bwci)
{
  unsigned long idx, addwidth = 0;

  for (idx = 0; idx < bwci->nextfreeBoundswithchar; idx++)
  {
    addwidth += bwci->spaceBoundswithchar[idx].rbound -
                bwci->spaceBoundswithchar[idx].lbound;
  }
  return addwidth;
}

void gt_fmindex_dfstraverse(const FMindex *fmindex,
                            unsigned int numofchars,
                            unsigned long totallength)
{
  GtArrayDfs_Boundsatdepth stack;
  GtArrayBoundswithchar bwci;
  Dfs_Boundsatdepth parent, child;
  unsigned long nonspecialwidth, parentwidth, idx, *rangeOccs;

  GT_INITARRAY(&stack,Dfs_Boundsatdepth);
  GT_INITARRAY(&bwci,Boundswithchar);
  child.lowerbound = 0;
  child.upperbound = totallength+1;
  child.depth = 0;
  child.endwithspecial = false;
  GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
  rangeOccs = gt_malloc(sizeof (*rangeOccs) * GT_MULT2(numofchars));
  while (stack.nextfreeDfs_Boundsatdepth > 0)
  {
    parent = stack.spaceDfs_Boundsatdepth[--stack.nextfreeDfs_Boundsatdepth];
    parentwidth = child.upperbound - parent.lowerbound;
    if (parent.endwithspecial)
    {
      printf("%lu special leaves\n",parentwidth);
    } else
    {
      if (parent.lowerbound < parent.upperbound)
      {
        gt_bwtrangesplitwithoutspecial(&bwci,rangeOccs,fmindex,
                                       parent.lowerbound,
                                       parent.upperbound);
        gt_assert(bwci.nextfreeBoundswithchar > 0);
        nonspecialwidth = dfsnonspecialwidth(&bwci);
        gt_assert(nonspecialwidth <= parentwidth);
        if (nonspecialwidth < parentwidth)
        {
          child.lowerbound = parent.lowerbound + nonspecialwidth;
          child.upperbound = parent.upperbound;
          child.depth = parent.depth + 1;
          child.endwithspecial = true;
          GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
        }
        gt_assert(bwci.spaceBoundswithchar != NULL);
        idx = bwci.nextfreeBoundswithchar - 1;
        while (true)
        {
          child.lowerbound = bwci.spaceBoundswithchar[idx].lbound;
          child.upperbound = bwci.spaceBoundswithchar[idx].rbound;
          child.depth = parent.depth + 1;
          child.endwithspecial = false;
          gt_assert(child.lowerbound < child.upperbound);
          GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
          if (idx > 0)
          {
            idx--;
          } else
          {
            break;
          }
        }
      } else
      {
        printf("leaf\nn");
      }
    }
  }
  GT_FREEARRAY(&stack,Dfs_Boundsatdepth);
  GT_FREEARRAY(&bwci,Boundswithchar);
  gt_free(rangeOccs);
}
