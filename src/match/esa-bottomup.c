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

#define TOP_ESA_BOTTOMUP   stackspace[nextfreeItvinfo-1]
#define POP_ESA_BOTTOMUP   stackspace[--nextfreeItvinfo]

#define PUSH_ESA_BOTTOMUP(LCP,LB,NOEDGE)\
        if (nextfreeItvinfo >= allocatedItvinfo)\
        {\
          gt_assert(nextfreeItvinfo == allocatedItvinfo);\
          stackspace = allocBUItvinfo(stackspace,\
                                      allocatedItvinfo,\
                                      allocatedItvinfo+incrementstacksize,\
                                      allocateBUinfo,\
                                      state);\
          allocatedItvinfo += incrementstacksize;\
        }\
        gt_assert(stackspace != NULL);\
        stackspace[nextfreeItvinfo].lcp = LCP;\
        stackspace[nextfreeItvinfo].lb = LB;\
        stackspace[nextfreeItvinfo].rb = GT_UNDEFBUBOUND;\
        stackspace[nextfreeItvinfo].noedge = NOEDGE;\
        nextfreeItvinfo++

#define SHOWBOOL_ESA_BOTTOMUP(B) ((B) ? '1' : '0')

typedef struct
{
  unsigned long lcp, lb, rb;
  bool noedge;
  GtBUinfo *info;
} GtBUItvinfo;

#ifdef SKDEBUG
static void showstack(const GtBUItvinfo *stackspace,
                      unsigned long nextfreeItvinfo)
{
  unsigned long idx;

  for (idx=0; idx<nextfreeItvinfo; idx++)
  {
    printf("# stack %lu: lcp=%lu,lb=%lu,noedge=%c\n",
            idx,
            stackspace[idx].lcp,
            stackspace[idx].lb,
            SHOWBOOL_ESA_BOTTOMUP(stackspace[idx].noedge));
  }
}
#endif

static GtBUItvinfo *allocBUItvinfo(GtBUItvinfo *ptr,
                                   unsigned long currentallocated,
                                   unsigned long allocated,
                                   GtBUinfo *(*allocateBUinfo)(GtBUstate *),
                                   GtBUstate *state)
{
  unsigned long idx;
  GtBUItvinfo *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof(*itvinfo) * allocated);
  gt_assert(allocated > currentallocated);
  for (idx=currentallocated; idx<allocated; idx++)
  {
    if (allocateBUinfo == NULL)
    {
      itvinfo[idx].info = NULL;
    } else
    {
      itvinfo[idx].info = allocateBUinfo(state);
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
    if (freeBUinfo != NULL)
    {
      freeBUinfo(ptr[idx].info,state);
    }
  }
  gt_free(ptr);
}

typedef struct
{
  unsigned long index,
                suffix;
} GtSuftabelem;

typedef struct
{
  GtSuftabelem leftelem,
               rightelem;
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
  gt_assert(suffix != GT_SUFTABELEMUNDEFFLAG);
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
    if (queuepair->leftelem.index < endpos)
    {
      printf("L %c %lu %lu %lu\n",SHOWBOOL_ESA_BOTTOMUP(firstedge),
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
    if (queuepair->rightelem.index < endpos)
    {
      printf("L %c %lu %lu %lu\n",SHOWBOOL_ESA_BOTTOMUP(firstedge),
             parent->lcp,
             parent->lb,
             queuepair->rightelem.suffix);
      GT_SUFTABELEMSETUNDEF(queuepair->rightelem);
      parent->noedge = false;
    }
  }
}

static void processbranchedge(bool firstedge,const GtBUItvinfo *fromitv,
                              const GtBUItvinfo *toitv)
{
  printf("B %c %lu %lu %lu %lu\n",SHOWBOOL_ESA_BOTTOMUP(firstedge),
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
  const unsigned long incrementstacksize = 32UL;
  GtQueuepair queuepair;
  GtBUItvinfo lastinterval, *stackspace = NULL;
  bool lastintervaldefined = false, firstedge, haserr = false;
  int retval;

  PUSH_ESA_BOTTOMUP(0,0,true);
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
    while (lcpvalue < TOP_ESA_BOTTOMUP.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP;
      lastinterval.rb = idx;
      lastintervaldefined = true;
      GtQueuepair_enumleaves(&lastinterval,&queuepair,idx+1);
      lb = lastinterval.lb;
      if (lcpvalue <= TOP_ESA_BOTTOMUP.lcp)
      {
        firstedge = TOP_ESA_BOTTOMUP.noedge;
        TOP_ESA_BOTTOMUP.noedge = false;
        processbranchedge(firstedge,&TOP_ESA_BOTTOMUP,&lastinterval);
        lastintervaldefined = false;
      }
    }
    if (lcpvalue > TOP_ESA_BOTTOMUP.lcp)
    {
      if (!lastintervaldefined)
      {
        GtQueuepair_enumleaves(&TOP_ESA_BOTTOMUP,&queuepair,lb);
        PUSH_ESA_BOTTOMUP(lcpvalue,lb,true);
      } else
      {
        PUSH_ESA_BOTTOMUP(lcpvalue,lb,false);
        processbranchedge(true,&TOP_ESA_BOTTOMUP,&lastinterval);
        lastintervaldefined = false;
      }
    } else
    {
      GtQueuepair_enumleaves(&TOP_ESA_BOTTOMUP,&queuepair,idx+1);
    }
  }
  lastinterval = POP_ESA_BOTTOMUP;
  lastinterval.rb = idx;
  lastintervaldefined = true;
  GtQueuepair_add(&queuepair,idx,idx);
  GtQueuepair_enumleaves(&lastinterval,&queuepair,idx+1);
  freeBUItvinfo(stackspace, allocatedItvinfo, freeBUinfo, state);
  return haserr ? -1 : 0;
}

int gt_runenumlcpvaluesBU(const char *inputindex,
                          GtLogger *logger,
                          GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(inputindex,
                                                   SARR_LCPTAB |
                                                   SARR_SUFTAB |
                                                   SARR_ESQTAB,
                                                   SEQ_scan,
                                                   err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr && gt_esa_bottomup(ssar, NULL, NULL, NULL, logger, err) != 0)
  {
    haserr = true;
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
