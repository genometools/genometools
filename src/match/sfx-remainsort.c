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

#include "core/chardef.h"
#include "core/queue.h"
#include "core/chardef.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "divmodmul.h"
#include "seqpos-def.h"
#include "sfx-remainsort.h"

typedef struct
{
  Seqpos key,
         suffixstart;
} Itventry;

typedef struct
{
  Seqpos *left,
         *right,
         depth;
} Pairsuffixptr;

struct Rmnsufinfo
{
  Seqpos *inversesuftab, *sortedsuffixes;
  unsigned long countovermaxdepthsingle;
  GtQueue *rangestobesorted;
  DefinedSeqpos previousdepth;
  Seqpos totallength;
};

Rmnsufinfo *initRmnsufinfo(Seqpos *sortedsuffixes,Seqpos totallength)
{
  Rmnsufinfo *rmnsufinfo;

  rmnsufinfo = gt_malloc(sizeof(Rmnsufinfo));
  rmnsufinfo->sortedsuffixes = sortedsuffixes;
  rmnsufinfo->countovermaxdepthsingle = 0;
  rmnsufinfo->rangestobesorted = gt_queue_new();
  rmnsufinfo->previousdepth.defined = false;
  rmnsufinfo->previousdepth.valueseqpos = 0;
  rmnsufinfo->totallength = totallength;
  return rmnsufinfo;
}

void addunsortedrange(Rmnsufinfo *rmnsufinfo,
                      Seqpos *left,Seqpos *right,Seqpos depth)
{
  if (left == right)
  {
    printf("already sorted(%lu,%lu,%lu)\n",
            (unsigned long) depth,
            (unsigned long) (left-rmnsufinfo->sortedsuffixes),
            (unsigned long) (right-rmnsufinfo->sortedsuffixes));
    rmnsufinfo->countovermaxdepthsingle++;
  } else
  {
    Pairsuffixptr *pairptr;

    pairptr = gt_malloc(sizeof(Pairsuffixptr));
    gt_assert(!rmnsufinfo->previousdepth.defined ||
              rmnsufinfo->previousdepth.valueseqpos <= depth);
    if (!rmnsufinfo->previousdepth.defined ||
        rmnsufinfo->previousdepth.valueseqpos < depth)
    {
      printf("new level with depth = %lu\n",(unsigned long) depth);
      rmnsufinfo->previousdepth.defined = true;
      rmnsufinfo->previousdepth.valueseqpos = depth;
    }
    pairptr->depth = depth;
    pairptr->left = left;
    pairptr->right = right;
    gt_queue_add(rmnsufinfo->rangestobesorted,pairptr);
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

static void setinversesuftab(Rmnsufinfo *rmnsufinfo,Seqpos idx,
                             Seqpos value)
{
  rmnsufinfo->inversesuftab[idx] = value;
}

static void inverserange(Rmnsufinfo *rmnsufinfo,Seqpos *left,Seqpos *right)
{
  Seqpos *ptr, startindex;

  startindex = (Seqpos) (left - rmnsufinfo->sortedsuffixes);
  for (ptr = left; ptr <= right; ptr++)
  {
    setinversesuftab(rmnsufinfo,*ptr,startindex);
  }
}

static void sortitv(Rmnsufinfo *rmnsufinfo,
                    Seqpos *left,Seqpos *right,Seqpos depth)
{
  Itventry *itvinfo;
  Seqpos startindex;
  unsigned long idx, rangestart;
  const unsigned long width = (unsigned long) (right - left + 1);

  itvinfo = gt_malloc(sizeof(Itventry) * width);
  for (idx=0; idx<width; idx++)
  {
    itvinfo[idx].suffixstart = left[idx];
    if (left[idx]+depth > rmnsufinfo->totallength)
    {
      fprintf(stderr,"left[%lu]+depth=%lu+%lu=%lu>%lu\n",
              idx,
              (unsigned long) left[idx],
              (unsigned long) depth,
              (unsigned long) (left[idx]+depth),
              (unsigned long) rmnsufinfo->totallength);
      exit(EXIT_FAILURE);
    }
    itvinfo[idx].key = rmnsufinfo->inversesuftab[left[idx]+depth];
  }
  qsort(itvinfo,(size_t) width,sizeof(Itventry),compareitv);
  for (idx=0; idx<width; idx++)
  {
    left[idx] = itvinfo[idx].suffixstart;
  }
  rangestart = 0;
  startindex = (Seqpos) (left - rmnsufinfo->sortedsuffixes);
  for (idx=1UL; idx<width; idx++)
  {
    if (itvinfo[idx-1].key != itvinfo[idx].key)
    {
      if (rangestart + 1 < idx)
      {
        addunsortedrange(rmnsufinfo,
                         left + rangestart,
                         left + idx - 1,
                         MULT2(depth));
        inverserange(rmnsufinfo,
                     left + rangestart,
                     left + idx - 1);
      } else
      {
        setinversesuftab(rmnsufinfo,left[rangestart],startindex+rangestart);
      }
      rangestart = idx;
    }
  }
  if (rangestart + 1 < width)
  {
    addunsortedrange(rmnsufinfo,
                     left + rangestart,
                     left + width - 1,
                     MULT2(depth));
    inverserange(rmnsufinfo,
                 left + rangestart,
                 left + width - 1);
  } else
  {
    setinversesuftab(rmnsufinfo,left[rangestart],startindex+rangestart);
  }
  gt_free(itvinfo);
}

static int putleftbound(void **elem,void *info, GT_UNUSED GtError *err)
{
  Pairsuffixptr *pairptr = *(Pairsuffixptr**) elem;

  inverserange((Rmnsufinfo *) info,pairptr->left,pairptr->right);
  return 0;
}

static void processRmnsufinfo(Rmnsufinfo *rmnsufinfo)
{
  Pairsuffixptr *pairptr;
  Seqpos idx;

  printf("# countovermaxdepth=%lu\n",
           gt_queue_size(rmnsufinfo->rangestobesorted));
  printf("# countovermaxdepthsingle=%lu\n",
          rmnsufinfo->countovermaxdepthsingle);
  rmnsufinfo->inversesuftab
    = gt_malloc(sizeof(Seqpos) * (rmnsufinfo->totallength+1));
  for (idx=0; idx <= rmnsufinfo->totallength; idx++)
  {
    rmnsufinfo->inversesuftab[rmnsufinfo->sortedsuffixes[idx]] = idx;
  }
  (void) gt_queue_iterate(rmnsufinfo->rangestobesorted,
                          putleftbound,
                          rmnsufinfo,
                          NULL);
  printf("# countovermaxdepth=%lu\n",
         gt_queue_size(rmnsufinfo->rangestobesorted));
  while (gt_queue_size(rmnsufinfo->rangestobesorted) > 0)
  {
    pairptr = gt_queue_get(rmnsufinfo->rangestobesorted);
    sortitv(rmnsufinfo,
            pairptr->left,
            pairptr->right,
            pairptr->depth);
    gt_free(pairptr);
  }
  gt_queue_delete(rmnsufinfo->rangestobesorted);
  rmnsufinfo->rangestobesorted = NULL;
  gt_free(rmnsufinfo->inversesuftab);
}

void wrapRmnsufinfo(Rmnsufinfo **rmnsufinfo)
{
  processRmnsufinfo(*rmnsufinfo);
  gt_free(*rmnsufinfo);
  *rmnsufinfo = NULL;
}
