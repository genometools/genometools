/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include "core/chardef.h"
#include "core/queue.h"
#include "core/chardef.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "divmodmul.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "sfx-remainsort.h"
#include "sfx-lcpsub.h"

typedef struct
{
  Seqpos key,
         suffixstart;
} Itventry;

typedef struct
{
  Seqpos *left,
         *right;
} Pairsuffixptr;

typedef struct
{
  Seqpos *left,
         *right,
         depth;
  unsigned long count, totalwidth;
} Firstwithnewdepth;

struct Rmnsufinfo
{
  Seqpos *inversesuftab,
         *presortedsuffixes;
  GtQueue *rangestobesorted;
  Seqpos partwidth,
         totallength,
         currentdepth;
  const Encodedsequence *encseq;
  Readmode readmode;
  unsigned long allocateditvinfo,
                currentqueuesize,
                maxqueuesize;
  Itventry *itvinfo;
  Firstwithnewdepth firstwithnewdepth;
  FILE *fplcptab, *fpllvtab;
};

Rmnsufinfo *newRmnsufinfo(Seqpos *presortedsuffixes,
                          const Encodedsequence *encseq,
                          Readmode readmode,
                          Seqpos partwidth,
                          FILE *fplcptab,
                          FILE *fpllvtab)
{
  Rmnsufinfo *rmnsufinfo;

  rmnsufinfo = gt_malloc(sizeof(Rmnsufinfo));
  rmnsufinfo->presortedsuffixes = presortedsuffixes;
  rmnsufinfo->rangestobesorted = gt_queue_new();
  rmnsufinfo->partwidth = partwidth;
  rmnsufinfo->totallength = getencseqtotallength(encseq);
  rmnsufinfo->encseq = encseq;
  rmnsufinfo->readmode = readmode;
  rmnsufinfo->allocateditvinfo = 0;
  rmnsufinfo->itvinfo = NULL;
  rmnsufinfo->currentqueuesize = 0;
  rmnsufinfo->maxqueuesize = 0;
  rmnsufinfo->firstwithnewdepth.depth = 0;
  rmnsufinfo->firstwithnewdepth.totalwidth = 0;
  rmnsufinfo->firstwithnewdepth.count = 0;
  rmnsufinfo->firstwithnewdepth.left = NULL;
  rmnsufinfo->firstwithnewdepth.right = NULL;
  rmnsufinfo->currentdepth = 0;
  rmnsufinfo->fplcptab = fplcptab;
  rmnsufinfo->fpllvtab = fpllvtab;
  return rmnsufinfo;
}

void addunsortedrange(Rmnsufinfo *rmnsufinfo,
                      Seqpos *left,Seqpos *right,Seqpos depth)
{
  Pairsuffixptr *pairptr;
  unsigned long width;

  gt_assert(left < right && depth > 0);
  gt_assert(rmnsufinfo->firstwithnewdepth.left == NULL ||
            (rmnsufinfo->firstwithnewdepth.depth > 0 &&
             rmnsufinfo->firstwithnewdepth.depth <= depth));
  if (rmnsufinfo->firstwithnewdepth.left == NULL ||
      rmnsufinfo->firstwithnewdepth.depth < depth)
  {
    if (rmnsufinfo->firstwithnewdepth.left != NULL)
    {
      printf("intervals in previous level=%lu (total=%lu,avg=%.2f,"
             "%.2f%% of all)\n",
              rmnsufinfo->firstwithnewdepth.count,
              rmnsufinfo->firstwithnewdepth.totalwidth,
              (double) rmnsufinfo->firstwithnewdepth.totalwidth/
                       rmnsufinfo->firstwithnewdepth.count,
              100.0 * (double) rmnsufinfo->firstwithnewdepth.totalwidth/
                               rmnsufinfo->totallength);
    }
    printf("enter new level with depth=%lu\n",(unsigned long) depth);
    rmnsufinfo->firstwithnewdepth.left = left;
    rmnsufinfo->firstwithnewdepth.right = right;
    rmnsufinfo->firstwithnewdepth.depth = depth;
    rmnsufinfo->firstwithnewdepth.count = 1UL;
    rmnsufinfo->firstwithnewdepth.totalwidth
      = (unsigned long) (right - left + 1);
  } else
  {
    rmnsufinfo->firstwithnewdepth.count++;
    rmnsufinfo->firstwithnewdepth.totalwidth
      += (unsigned long) (right - left + 1);
  }
  width = (unsigned long) (right - left + 1);
  if (rmnsufinfo->allocateditvinfo < width)
  {
    gt_assert(rmnsufinfo->itvinfo == NULL);
    rmnsufinfo->allocateditvinfo = width;
  }
  pairptr = gt_malloc(sizeof(Pairsuffixptr));
  pairptr->left = left;
  pairptr->right = right;
  gt_queue_add(rmnsufinfo->rangestobesorted,pairptr);
  rmnsufinfo->currentqueuesize++;
  if (rmnsufinfo->maxqueuesize < rmnsufinfo->currentqueuesize)
  {
    rmnsufinfo->maxqueuesize = rmnsufinfo->currentqueuesize;
  }
}

static int compareitv(const void *a,const void *b)
{
  const Itventry *itva = (const Itventry *) a,
                 *itvb = (const Itventry *) b;

  if (itva->key < itvb->key)
  {
    return -1;
  }
  if (itva->key > itvb->key)
  {
    return 1;
  }
  return 0;
}

static void inverserange(Rmnsufinfo *rmnsufinfo,Seqpos *left,Seqpos *right)
{
  Seqpos *ptr, startindex;

  startindex = (Seqpos) (left - rmnsufinfo->presortedsuffixes);
  for (ptr = left; ptr <= right; ptr++)
  {
    rmnsufinfo->inversesuftab[*ptr] = startindex;
  }
}

static void sortitv(Rmnsufinfo *rmnsufinfo,Seqpos *left,Seqpos *right)
{
  Seqpos startindex;
  unsigned long idx, rangestart;
  const unsigned long width = (unsigned long) (right - left + 1);

  if (rmnsufinfo->firstwithnewdepth.left == left &&
      rmnsufinfo->firstwithnewdepth.right == right)
  {
    rmnsufinfo->currentdepth = rmnsufinfo->firstwithnewdepth.depth;
  }
  gt_assert(rmnsufinfo->allocateditvinfo >= width);
  for (idx=0; idx<width; idx++)
  {
    rmnsufinfo->itvinfo[idx].suffixstart = left[idx];
    if (left[idx]+rmnsufinfo->currentdepth > rmnsufinfo->totallength)
    {
      fprintf(stderr,"left[%lu]+depth=%lu+%lu=%lu>%lu\n",
              idx,
              (unsigned long) left[idx],
              (unsigned long) rmnsufinfo->currentdepth,
              (unsigned long) (left[idx]+rmnsufinfo->currentdepth),
              (unsigned long) rmnsufinfo->totallength);
      exit(EXIT_FAILURE); /* programm error */
    }
    rmnsufinfo->itvinfo[idx].key
      = rmnsufinfo->inversesuftab[left[idx]+rmnsufinfo->currentdepth];
  }
  qsort(rmnsufinfo->itvinfo,(size_t) width,sizeof(Itventry),compareitv);
  for (idx=0; idx<width; idx++)
  {
    left[idx] = rmnsufinfo->itvinfo[idx].suffixstart;
  }
  rangestart = 0;
  startindex = (Seqpos) (left - rmnsufinfo->presortedsuffixes);
  for (idx=1UL; idx<width; idx++)
  {
    if (rmnsufinfo->itvinfo[idx-1].key != rmnsufinfo->itvinfo[idx].key)
    {
      if (rangestart + 1 < idx)
      {
        addunsortedrange(rmnsufinfo,
                         left + rangestart,
                         left + idx - 1,
                         MULT2(rmnsufinfo->currentdepth));
        inverserange(rmnsufinfo,
                     left + rangestart,
                     left + idx - 1);
      } else
      {
        rmnsufinfo->inversesuftab[left[rangestart]] = startindex+rangestart;
      }
      rangestart = idx;
    }
  }
  if (rangestart + 1 < width)
  {
    addunsortedrange(rmnsufinfo,
                     left + rangestart,
                     left + width - 1,
                     MULT2(rmnsufinfo->currentdepth));
    inverserange(rmnsufinfo,
                 left + rangestart,
                 left + width - 1);
  } else
  {
    rmnsufinfo->inversesuftab[left[rangestart]] = startindex+rangestart;
  }
}

static int putleftbound(void **elem,void *info, GT_UNUSED GtError *err)
{
  Pairsuffixptr *pairptr = *(Pairsuffixptr**) elem;

  inverserange((Rmnsufinfo *) info,pairptr->left,pairptr->right);
  return 0;
}

static void sortremainingsuffixes(Rmnsufinfo *rmnsufinfo)
{
  Pairsuffixptr *pairptr;
  Seqpos idx;

  printf("# countovermaxdepth=%lu\n",
           gt_queue_size(rmnsufinfo->rangestobesorted));
  rmnsufinfo->inversesuftab
    = gt_malloc(sizeof(Seqpos) * (rmnsufinfo->totallength+1));
  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    rmnsufinfo->inversesuftab[rmnsufinfo->presortedsuffixes[idx]] = idx;
  }
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialidx;

    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    specialidx = rmnsufinfo->partwidth;
    while (nextspecialrangeiterator(&range,sri))
    {
      for (idx = range.leftpos; idx < range.rightpos; idx++)
      {
        rmnsufinfo->inversesuftab[idx] = specialidx++;
      }
    }
    gt_assert(specialidx == rmnsufinfo->totallength);
    freespecialrangeiterator(&sri);
  }
  /* now we can unmap the encoded sequence */
  rmnsufinfo->inversesuftab[rmnsufinfo->totallength] = rmnsufinfo->totallength;
  (void) gt_queue_iterate(rmnsufinfo->rangestobesorted,
                          putleftbound,
                          rmnsufinfo,
                          NULL);
  rmnsufinfo->itvinfo
    = gt_malloc(sizeof(Itventry) * (rmnsufinfo->allocateditvinfo));
  gt_assert(rmnsufinfo->currentqueuesize  ==
            gt_queue_size(rmnsufinfo->rangestobesorted));
  while (gt_queue_size(rmnsufinfo->rangestobesorted) > 0)
  {
    pairptr = gt_queue_get(rmnsufinfo->rangestobesorted);
    gt_assert(rmnsufinfo->currentqueuesize > 0);
    rmnsufinfo->currentqueuesize--;
    sortitv(rmnsufinfo,
            pairptr->left,
            pairptr->right);
    gt_free(pairptr);
  }
  printf("maxqueuesize = %lu\n",rmnsufinfo->maxqueuesize);
}

static void lineartimelcpcomputation(const Rmnsufinfo *rmnsufinfo)
{
  Seqpos idx, h = 0, *lcptab;

#ifndef NDEBUG
  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    gt_assert(rmnsufinfo->inversesuftab[rmnsufinfo->presortedsuffixes[idx]]
              == idx);
  }
#endif
  lcptab = gt_malloc(sizeof(Seqpos) * rmnsufinfo->partwidth);
  lcptab[0] = 0;
  for (idx=0; idx <= rmnsufinfo->totallength; idx++)
  {
    Seqpos k = rmnsufinfo->inversesuftab[idx];
    if (k > 0 && k < rmnsufinfo->partwidth)
    {
      Seqpos jdx = rmnsufinfo->presortedsuffixes[k-1];
      while (idx+h < rmnsufinfo->totallength &&
             jdx+h < rmnsufinfo->totallength)
      {
        Uchar cc1, cc2;

        cc1 = getencodedchar(rmnsufinfo->encseq,idx+h,rmnsufinfo->readmode);
        cc2 = getencodedchar(rmnsufinfo->encseq,jdx+h,rmnsufinfo->readmode);
        if (cc1 == cc2 && ISNOTSPECIAL(cc1))
        {
          h++;
        } else
        {
          break;
        }
      }
      lcptab[k] = h;
    }
    if (h > 0)
    {
      h--;
    }
  }
  multioutlcpvalues(lcptab,(unsigned long) rmnsufinfo->partwidth,
                    rmnsufinfo->fplcptab,rmnsufinfo->fpllvtab);
  gt_free(lcptab);
}

void wrapRmnsufinfo(Rmnsufinfo **rmnsufinfoptr)
{
  Rmnsufinfo *rmnsufinfo = *rmnsufinfoptr;

  sortremainingsuffixes(rmnsufinfo);
  if (rmnsufinfo->fplcptab != NULL)
  {
    gt_assert(rmnsufinfo->fpllvtab != NULL);
    lineartimelcpcomputation(rmnsufinfo);
  }
  gt_free(rmnsufinfo->itvinfo);
  rmnsufinfo->itvinfo = NULL;
  gt_queue_delete(rmnsufinfo->rangestobesorted);
  rmnsufinfo->rangestobesorted = NULL;
  gt_free(rmnsufinfo->inversesuftab);
  rmnsufinfo->inversesuftab = NULL;
  gt_free(rmnsufinfo);
  rmnsufinfoptr = NULL;
}
