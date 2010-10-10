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
  union
  {
    struct
    {
      unsigned long lowerbound,
                    upperbound;
    } pckitv;
    unsigned long remainingspecial;
  } either;
  unsigned long depth;
  bool isinterval;
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

#ifdef SKDEBUG
static void showDfs_Boundsatdepth(const Dfs_Boundsatdepth *bd)
{
  if (bd->isinterval)
  {
    printf("pop l=%lu u=%lu d=%lu\n",bd->either.pckitv.lowerbound,
                                     bd->either.pckitv.upperbound,
                                     bd->depth);
  } else
  {
    printf("pop w=%lu d=%lu\n",bd->either.remainingspecial,bd->depth);
  }
}
#endif

void gt_fmindex_dfstraverse(const FMindex *fmindex,
                            unsigned int numofchars,
                            unsigned long totallength)
{
  GtArrayDfs_Boundsatdepth stack;
  GtArrayBoundswithchar bwci;
  Boundswithchar *bwciptr;
  Dfs_Boundsatdepth parent, child;
  unsigned long nonspecialwidth, parentwidth, *rangeOccs;
#undef OUTPUT
#ifdef OUTPUT
  bool firstleaf = true;
#endif

  GT_INITARRAY(&stack,Dfs_Boundsatdepth);
  bwci.spaceBoundswithchar = gt_malloc(sizeof (*bwci.spaceBoundswithchar) *
                                       (numofchars+1));
  bwci.nextfreeBoundswithchar = 0;
  bwci.allocatedBoundswithchar = (unsigned long) (numofchars+1);
  child.isinterval = true;
  child.depth = 0;
  child.either.pckitv.lowerbound = 0;
  child.either.pckitv.upperbound = totallength+1;
  GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
  rangeOccs = gt_malloc(sizeof (*rangeOccs) * GT_MULT2(numofchars));
  while (stack.nextfreeDfs_Boundsatdepth > 0)
  {
    parent = stack.spaceDfs_Boundsatdepth[--stack.nextfreeDfs_Boundsatdepth];
#ifdef SKDEBUG
    showDfs_Boundsatdepth(&parent);
#endif
    if (!parent.isinterval)
    {
#ifdef OUTPUT
      printf("%lu special leaves with lcp %lu\n",
               parent.either.remainingspecial,
               parent.depth);
#endif
    } else
    {
      gt_assert(parent.either.pckitv.lowerbound <
                parent.either.pckitv.upperbound);
      parentwidth = parent.either.pckitv.upperbound -
                    parent.either.pckitv.lowerbound;
      if (parentwidth == 1UL)
      {
#ifdef OUTPUT
        printf("leaf");
        if (!firstleaf)
        {
          printf(" with lcp %lu\n",parent.depth);
        } else
        {
          printf("\n");
          firstleaf = false;
        }
#endif
      } else
      {
        gt_assert(parentwidth >= 2UL);
        gt_bwtrangesplitwithoutspecial(&bwci,rangeOccs,fmindex,
                                       parent.either.pckitv.lowerbound,
                                       parent.either.pckitv.upperbound);
        if (bwci.nextfreeBoundswithchar == 1UL && 
            parentwidth == bwci.spaceBoundswithchar[0].rbound - 
                           bwci.spaceBoundswithchar[0].lbound)
        {
          child.isinterval = true;
          child.depth = parent.depth + 1;
          child.either.pckitv.lowerbound = bwci.spaceBoundswithchar[0].lbound;
          child.either.pckitv.upperbound = bwci.spaceBoundswithchar[0].rbound;
          gt_assert(child.either.pckitv.lowerbound <
                    child.either.pckitv.upperbound);
          GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
        } else
        {
          nonspecialwidth = dfsnonspecialwidth(&bwci);
          gt_assert(nonspecialwidth <= parentwidth);
          if (nonspecialwidth < parentwidth)
          {
            child.isinterval = false;
            child.depth = parent.depth + 1;
            child.either.remainingspecial = parentwidth - nonspecialwidth;
            GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
          }
          gt_assert(bwci.spaceBoundswithchar != NULL);
          for (bwciptr = bwci.spaceBoundswithchar+bwci.nextfreeBoundswithchar-1;
               bwciptr >= bwci.spaceBoundswithchar;
               bwciptr--)
          {
            child.isinterval = true;
            child.depth = parent.depth + 1;
            child.either.pckitv.lowerbound = bwciptr->lbound;
            child.either.pckitv.upperbound = bwciptr->rbound;
            gt_assert(child.either.pckitv.lowerbound <
                      child.either.pckitv.upperbound);
            GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
          }
        }
      }
    }
  }
  GT_FREEARRAY(&stack,Dfs_Boundsatdepth);
  GT_FREEARRAY(&bwci,Boundswithchar);
  gt_free(rangeOccs);
}
