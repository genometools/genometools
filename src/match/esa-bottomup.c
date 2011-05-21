/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/unused_api.h"
#include "core/ma.h"
#include "esa-bottomup.h"
#include "esa-seqread.h"

#define TOP           stackspace[nextfreeItvinfo-1]
#define INCSTACKSIZE  32

#define PUSHBOTTOMUP(LCP,LB,NOEDGE)\
        if (nextfreeItvinfo >= allocatedItvinfo)\
        {\
          gt_assert(nextfreeItvinfo == allocatedItvinfo);\
          stackspace = allocBUItvinfo(stackspace,\
                                      allocatedItvinfo,\
                                      allocatedItvinfo+INCSTACKSIZE,\
                                      allocateBUinfo,\
                                      state);\
          allocatedItvinfo += INCSTACKSIZE;\
        }\
        gt_assert(stackspace != NULL);\
        stackspace[nextfreeItvinfo].lcp = LCP;\
        stackspace[nextfreeItvinfo].lb = LB;\
        stackspace[nextfreeItvinfo].rb = GT_UNDEFBUBOUND;\
        stackspace[nextfreeItvinfo].noedge = NOEDGE;\
        nextfreeItvinfo++

typedef struct
{
  unsigned long lcp, lb, rb;
  bool noedge;
  GtBUinfo *info;
} GtBUItvinfo;

static GtBUItvinfo *allocBUItvinfo(GtBUItvinfo *ptr,
                                   unsigned long currentallocated,
                                   unsigned long allocated,
                                   GtBUinfo *(*allocateBUinfo)(GtBUstate *),
                                   GtBUstate *state)
{
  unsigned long i;
  GtBUItvinfo *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof(*itvinfo) * allocated);
  if (allocateBUinfo != NULL)
  {
    gt_assert(allocated > currentallocated);
    for (i=currentallocated; i<allocated; i++)
    {
      itvinfo[i].info = allocateBUinfo(state);
    }
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}

static void freeBUItvinfo(GtBUItvinfo *ptr,
                          unsigned long allocated,
                          void (*freeBUinfo)(GtBUinfo *,GtBUstate *),
                          GtBUstate *state)
{
  unsigned long idx;

  for (idx=0; idx<allocated; idx++)
  {
    freeBUinfo(ptr[idx].info,state);
  }
  gt_free(ptr);
}

typedef struct
{
  unsigned long index, suffix;
} GtSuftabelem;

typedef struct
{
  GtSuftabelem leftelem, rightelem;
} GtQueuepair;

#define GT_SUFTABELEMUNDEFFLAG  ULONG_MAX
#define GT_SUFTABELEMISUNDEF(X) ((X).suffix == GT_SUFTABELEMUNDEFFLAG ? true\
                                                                      : false)
#define GT_SUFTABELEMSETUNDEF(X) (X).suffix = GT_SUFTABELEMUNDEFFLAG
#define GT_UNDEFBUBOUND          ULONG_MAX

static void GtQueuepair_init(GtQueuepair *queuepair)
{
  queuepair->leftelem.index = queuepair->rightelem.index = 0;
  GT_SUFTABELEMSETUNDEF(queuepair->leftelem);
  GT_SUFTABELEMSETUNDEF(queuepair->rightelem);
}

static void GtQueuepair_add(GtQueuepair *queuepair,unsigned long index,
                            unsigned long suffix)
{
  gt_assert(GT_SUFTABELEMISUNDEF(queuepair->leftelem));
  queuepair->leftelem = queuepair->rightelem;
  queuepair->rightelem.index = index;
  queuepair->rightelem.suffix = suffix;
}

static void GtQueuepair_enumleaves(GtBUItvinfo *parent,
                                   GtQueuepair *queuepair,
                                   unsigned long endpos)
{
  bool firstedge;

  gt_assert(parent != NULL);
  firstedge = parent->noedge;
  if (!GT_SUFTABELEMISUNDEF(queuepair->leftelem))
  {
    if (queuepair->leftelem.index <= endpos)
    {
      printf("L %s %lu %lu %lu\n",firstedge ? "true" : "false",
             parent->lcp,parent->lb,queuepair->leftelem.suffix);
      GT_SUFTABELEMSETUNDEF(queuepair->leftelem);
      firstedge = false;
      parent->noedge = false;
    } else
    {
      return;
    }
  }
  if (!GT_SUFTABELEMISUNDEF(queuepair->rightelem))
  {
    if (queuepair->rightelem.index <= endpos)
    {
      printf("L %s %lu %lu %lu\n",firstedge ? "true" : "false",
             parent->lcp,parent->lb,queuepair->rightelem.suffix);
      GT_SUFTABELEMSETUNDEF(queuepair->rightelem);
    }
  }
}

static void processbranchedge(bool firstedge,const GtBUItvinfo *fromitv,
                              const GtBUItvinfo *toitv)
{
  printf("B %c %lu %lu %lu %lu\n",firstedge ? '1' : '0',
         fromitv->lcp,fromitv->lb,toitv->lcp,toitv->lb);
}

int gt_esa_bottomup(Sequentialsuffixarrayreader *ssar,
                    GtBUinfo *(*allocateBUinfo)(GtBUstate *),
                    void(*freeBUinfo)(GtBUinfo *,GtBUstate *),
                    GtBUstate *state,
                    GT_UNUSED GtLogger *logger,
                    GtError *err)
{
  unsigned long lcpvalue,
                previoussuffix,
                lb,
                idx,
                allocatedItvinfo = 0,
                nextfreeItvinfo = 0;
  GtBUItvinfo lastinterval, *stackspace = NULL;
  GtQueuepair queuepair;
  bool firstedge, lastintervaldefined = false, haserr = false;
  int retval;

  PUSHBOTTOMUP(0,0,true);
  GtQueuepair_init(&queuepair);
  for (idx = 0; !haserr; idx++)
  {
    retval = gt_nextSequentiallcpvalue(&lcpvalue,ssar,err);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    retval = gt_nextSequentialsuftabvalue(&previoussuffix,ssar);
    gt_assert(retval >= 0);
    if (retval == 0)
    {
      gt_error_set(err,"Missing value in suftab");
      haserr = true;
      break;
    }
    lb = idx;
    GtQueuepair_add(&queuepair,idx,previoussuffix);
    while (lcpvalue < TOP.lcp)
    {
      lastinterval = stackspace[nextfreeItvinfo--];
      lastinterval.rb = idx;
      lastintervaldefined = true;
      GtQueuepair_enumleaves(&lastinterval,&queuepair,idx);
      lb = lastinterval.lb;
      if (lcpvalue <= TOP.lcp)
      {
        firstedge = TOP.noedge;
        TOP.noedge = false;
        processbranchedge(firstedge,&TOP,&lastinterval);
        lastintervaldefined = false;
      }
    }
    if (lcpvalue > TOP.lcp)
    {
      if (!lastintervaldefined)
      {
        GtQueuepair_enumleaves(&TOP,&queuepair,lb-1);
        PUSHBOTTOMUP(lcpvalue,lb,true);
      } else
      {
        PUSHBOTTOMUP(lcpvalue,lb,false);
        processbranchedge(true,&TOP,&lastinterval);
        lastintervaldefined = false;
      }
    } else
    {
      GtQueuepair_enumleaves(&TOP,&queuepair,idx);
    }
  }
  lastinterval = stackspace[nextfreeItvinfo--];
  lastinterval.rb = idx;
  lastintervaldefined = true;
  GtQueuepair_add(&queuepair,idx,idx);
  GtQueuepair_enumleaves(&lastinterval,&queuepair,idx);
  freeBUItvinfo(stackspace, allocatedItvinfo, freeBUinfo, state);
  return haserr ? -1 : 0;
}
