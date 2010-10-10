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
      GtUchar inchar;
    } pckitv;
    unsigned long remainingspecial;
  } either;
  unsigned long depth,
                lcpnodedepth;
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

#undef SKDEBUG
#ifdef SKDEBUG
static void showcurrentpath(const GtArrayGtUchar *currentpath)
{
  unsigned long idx;

  printf("path=");
  for (idx = 0; idx < currentpath->nextfreeGtUchar; idx++)
  {
    printf("%u ",(unsigned int) currentpath->spaceGtUchar[idx]);
  }
  printf("\n");
}

static void showDfs_Boundsatdepth(const char *kind,const Dfs_Boundsatdepth *bd)
{
  printf("%s ",kind);
  if (bd->isinterval)
  {
    printf("l=%lu u=%lu i=%u d=%lu\n",bd->either.pckitv.lowerbound,
                                      bd->either.pckitv.upperbound,
                                      (unsigned int)
                                      bd->either.pckitv.inchar,
                                      bd->depth);
  } else
  {
    printf("w=%lu d=%lu\n",bd->either.remainingspecial,bd->depth);
  }
}
#endif

int gt_fmindex_dfstraverse(const FMindex *fmindex,
                           unsigned int numofchars,
                           unsigned long totallength,
                           Processlcp processlcp,
                           void *processlcpdata,
                           GtError *err)
{
  GtArrayDfs_Boundsatdepth stack;
  GtArrayBoundswithchar bwci;
  Boundswithchar *bwciptr;
  Dfs_Boundsatdepth parent, child;
  GtArrayGtUchar currentpath;
  unsigned long nonspecialwidth, parentwidth, *rangeOccs;
  bool haserr = false, firstleaf = true;

  GT_INITARRAY(&stack,Dfs_Boundsatdepth);
  GT_INITARRAY(&currentpath,GtUchar);
  bwci.spaceBoundswithchar = gt_malloc(sizeof (*bwci.spaceBoundswithchar) *
                                       (numofchars+1));
  bwci.nextfreeBoundswithchar = 0;
  bwci.allocatedBoundswithchar = (unsigned long) (numofchars+1);
  child.isinterval = true;
  child.depth = 0;
  child.lcpnodedepth = 0;
  child.either.pckitv.lowerbound = 0;
  child.either.pckitv.upperbound = totallength+1;
  child.either.pckitv.inchar = (GtUchar) numofchars; /* undefined */
#ifdef SKDEBUG
  showDfs_Boundsatdepth("push",&child);
#endif
  GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
  rangeOccs = gt_malloc(sizeof (*rangeOccs) * GT_MULT2(numofchars));
  while (!haserr && stack.nextfreeDfs_Boundsatdepth > 0)
  {
    parent = stack.spaceDfs_Boundsatdepth[--stack.nextfreeDfs_Boundsatdepth];
#ifdef SKDEBUG
    showDfs_Boundsatdepth("pop",&parent);
#endif
    if (!parent.isinterval)
    {

      gt_assert(parent.depth > 0);
      if (processlcp != NULL)
      {
        unsigned long idx;

        for (idx = 0; idx<parent.either.remainingspecial; idx++)
        {
          if (!firstleaf)
          {
            if (processlcp(processlcpdata,(idx == 0)
                                          ? parent.lcpnodedepth
                                          : parent.depth - 1,
                                          err) != 0)
            {
              haserr = true;
              break;
            }
          } else
          {
            firstleaf = false;
          }
        }
      }
    } else
    {
      gt_assert(parent.either.pckitv.lowerbound <
                parent.either.pckitv.upperbound);
      if (parent.depth > 0)
      {
        gt_assert(parent.either.pckitv.inchar < (GtUchar) numofchars);
        if (parent.depth - 1 >= currentpath.allocatedGtUchar)
        {
          currentpath.allocatedGtUchar += 32UL;
          currentpath.spaceGtUchar
            = gt_realloc(currentpath.spaceGtUchar,
                         sizeof (*currentpath.spaceGtUchar) *
                         currentpath.allocatedGtUchar);
        }
        gt_assert(currentpath.spaceGtUchar != NULL);
        currentpath.spaceGtUchar[parent.depth - 1]
          = parent.either.pckitv.inchar;
        currentpath.nextfreeGtUchar = parent.depth;
#ifdef SKDEBUG
        showcurrentpath(&currentpath);
#endif
      }
      parentwidth = parent.either.pckitv.upperbound -
                    parent.either.pckitv.lowerbound;
      if (parentwidth == 1UL)
      {
        if (!firstleaf)
        {
          if (processlcp != NULL)
          {
            if (processlcp(processlcpdata,parent.lcpnodedepth,err) != 0)
            {
              haserr = true;
            }
          }
        } else
        {
          firstleaf = false;
        }
      } else
      {
        gt_assert(parentwidth >= 2UL);
        gt_bwtrangesplitwithoutspecial(&bwci,rangeOccs,fmindex,
                                       parent.either.pckitv.lowerbound,
                                       parent.either.pckitv.upperbound);
        nonspecialwidth = dfsnonspecialwidth(&bwci);
#ifdef SKDEBUG
        printf("split %lu %lu into %lu intervals of width %lu\n",
               parent.either.pckitv.lowerbound,
               parent.either.pckitv.upperbound,
               bwci.nextfreeBoundswithchar,
               nonspecialwidth);
#endif
        gt_assert(nonspecialwidth <= parentwidth);
        if (nonspecialwidth < parentwidth)
        {
          child.isinterval = false;
          child.depth = parent.depth + 1;
          child.either.remainingspecial = parentwidth - nonspecialwidth;
          if (bwci.nextfreeBoundswithchar > 0)
          {
            child.lcpnodedepth = parent.depth;
          } else
          {
            child.lcpnodedepth = parent.lcpnodedepth;
          }
#ifdef SKDEBUG
          showDfs_Boundsatdepth("special push",&child);
#endif
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
          if (bwciptr > bwci.spaceBoundswithchar)
          {
            child.lcpnodedepth = parent.depth;
          } else
          {
            child.lcpnodedepth = parent.lcpnodedepth;
          }
          gt_assert(bwciptr->inchar < (GtUchar) numofchars);
          child.either.pckitv.inchar = bwciptr->inchar;
          gt_assert(child.either.pckitv.lowerbound <
                    child.either.pckitv.upperbound);
#ifdef SKDEBUG
          showDfs_Boundsatdepth("push",&child);
#endif
          GT_STOREINARRAY(&stack,Dfs_Boundsatdepth,128,child);
        }
      }
    }
  }
  GT_FREEARRAY(&stack,Dfs_Boundsatdepth);
  GT_FREEARRAY(&bwci,Boundswithchar);
  GT_FREEARRAY(&currentpath,GtUchar);
  gt_free(rangeOccs);
  return haserr ? -1 : 0;
}
