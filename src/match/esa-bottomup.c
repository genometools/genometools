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
#include "core/ma.h"
#include "esa-bottomup.h"
#include "esa-seqread.h"

#define TOP_ESA_BOTTOMUP   stackspace[nextfreeItvinfo-1]
#define POP_ESA_BOTTOMUP   (stackspace + (--nextfreeItvinfo))

#define PUSH_ESA_BOTTOMUP(LCP,LB)\
        if (nextfreeItvinfo >= allocatedItvinfo)\
        {\
          gt_assert(nextfreeItvinfo == allocatedItvinfo);\
          stackspace = allocateBUstack(stackspace,\
                                       allocatedItvinfo,\
                                       allocatedItvinfo+incrementstacksize,\
                                       allocateBUinfo,\
                                       bustate);\
          allocatedItvinfo += incrementstacksize;\
        }\
        gt_assert(stackspace != NULL);\
        stackspace[nextfreeItvinfo].lcp = LCP;\
        stackspace[nextfreeItvinfo].lb = LB;\
        stackspace[nextfreeItvinfo].rb = ULONG_MAX;\
        nextfreeItvinfo++

typedef struct
{
  unsigned long lcp, lb, rb;
  GtBUinfo *info;
} GtBUItvinfo;

#ifdef SKDEBUG
static void showstack(const GtBUItvinfo *stackspace,
                      unsigned long nextfreeItvinfo)
{
  unsigned long idx;

  for (idx=0; idx<nextfreeItvinfo; idx++)
  {
    printf("# stack %lu: lcp=%lu,lb=%lu\n",
            idx,
            stackspace[idx].lcp,
            stackspace[idx].lb);
  }
}
#endif

static GtBUItvinfo *allocateBUstack(GtBUItvinfo *ptr,
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

int gt_esa_bottomup(Sequentialsuffixarrayreader *ssar,
                    GtBUinfo *(*allocateBUinfo)(GtBUstate *),
                    void(*freeBUinfo)(GtBUinfo *,GtBUstate *),
                    int (*processleafedge)(bool,
                                           unsigned long,
                                           unsigned long,
                                           GtBUinfo *,
                                           unsigned long,
                                           GtBUstate *,
                                           GtError *err),
                    int (*processbranchingedge)(bool firstsucc,
                                                unsigned long,
                                                unsigned long,
                                                GtBUinfo *,
                                                unsigned long,
                                                unsigned long,
                                                unsigned long,
                                                GtBUinfo *,
                                                GtBUstate *,
                                                GtError *),
                    GtBUstate *bustate,
                    GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix,
                idx,
                nonspecials,
                allocatedItvinfo = 0,
                nextfreeItvinfo = 0;
  GtBUItvinfo *lastinterval = NULL, *stackspace = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;
  int retval;

  PUSH_ESA_BOTTOMUP(0,0);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  for (idx = 0; idx < nonspecials; idx++)
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
    if (lcpvalue <= TOP_ESA_BOTTOMUP.lcp)
    {
      if (TOP_ESA_BOTTOMUP.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge(firstedge,
                          TOP_ESA_BOTTOMUP.lcp,
                          TOP_ESA_BOTTOMUP.lb,
                          TOP_ESA_BOTTOMUP.info,
                          previoussuffix,bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP;
      lastinterval->rb = idx;
      if (lcpvalue <= TOP_ESA_BOTTOMUP.lcp)
      {
        gt_assert(lastinterval->info == NULL ||
                  lastinterval->info != TOP_ESA_BOTTOMUP.info);
        if (TOP_ESA_BOTTOMUP.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
        if (processbranchingedge(firstedge,
                                 TOP_ESA_BOTTOMUP.lcp,
                                 TOP_ESA_BOTTOMUP.lb,
                                 TOP_ESA_BOTTOMUP.info,
                                 lastinterval->lcp,
                                 lastinterval->lb,
                                 lastinterval->rb,
                                 lastinterval->info,
                                 bustate,
                                 err) != 0)
        {
          haserr = true;
          break;
        }
        lastinterval = NULL;
      }
    }
    if (haserr)
    {
      break;
    }
    if (lcpvalue > TOP_ESA_BOTTOMUP.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervallb = lastinterval->lb,
                      lastintervalrb = lastinterval->rb;
        PUSH_ESA_BOTTOMUP(lcpvalue,lastintervallb);
        if (processbranchingedge(true,
                                 TOP_ESA_BOTTOMUP.lcp,
                                 TOP_ESA_BOTTOMUP.lb,
                                 TOP_ESA_BOTTOMUP.info,
                                 lastintervallcp,
                                 lastintervallb,
                                 lastintervalrb,
                                 NULL,
                                 bustate,err) != 0)
        {
          haserr = true;
          break;
        }
        lastinterval = NULL;
      } else
      {
        PUSH_ESA_BOTTOMUP(lcpvalue,idx);
        if (processleafedge(true,
                            TOP_ESA_BOTTOMUP.lcp,
                            TOP_ESA_BOTTOMUP.lb,
                            TOP_ESA_BOTTOMUP.info,
                            previoussuffix,bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  freeBUItvinfo(stackspace, allocatedItvinfo, freeBUinfo, bustate);
  return haserr ? -1 : 0;
}
