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
#include "esa-seqread.h"

static void initBUinfo_rdjce(BUinfo_rdjce *,BUstate_rdjce *);

static void freeBUinfo_rdjce(BUinfo_rdjce *,BUstate_rdjce *);

static int processleafedge_rdjce(bool,
    unsigned long,
    unsigned long,
    BUinfo_rdjce *,
    unsigned long,
    BUstate_rdjce *,
    GtError *err);

static int processbranchingedge_rdjce(bool firstsucc,
    unsigned long,
    unsigned long,
    BUinfo_rdjce *,
    unsigned long,
    unsigned long,
    unsigned long,
    BUinfo_rdjce *,
    BUstate_rdjce *,
    GtError *);

static int processlcpinterval_rdjce(unsigned long,
    unsigned long,
    unsigned long,
    BUinfo_rdjce *,
    BUstate_rdjce *,
    GtError *err);

#define TOP_ESA_BOTTOMUP_rdjce   stackspace[nextfreeItvinfo-1]
#define POP_ESA_BOTTOMUP_rdjce   (stackspace + (--nextfreeItvinfo))

#define PUSH_ESA_BOTTOMUP_rdjce(LCP,LB)\
        if (nextfreeItvinfo >= allocatedItvinfo)\
        {\
          gt_assert(nextfreeItvinfo == allocatedItvinfo);\
          stackspace = allocateBUstack_rdjce(stackspace,\
                                       allocatedItvinfo,\
                                       allocatedItvinfo+incrementstacksize,\
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
  BUinfo_rdjce info;
} GtBUItvinfo_rdjce;

static GtBUItvinfo_rdjce *allocateBUstack_rdjce(GtBUItvinfo_rdjce *ptr,
                                   unsigned long currentallocated,
                                   unsigned long allocated,
                                   BUstate_rdjce *state)
{
  unsigned long idx;
  GtBUItvinfo_rdjce *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof(*itvinfo) * allocated);
  gt_assert(allocated > currentallocated);
  for (idx=currentallocated; idx<allocated; idx++)
  {
    initBUinfo_rdjce(&itvinfo[idx].info,state);
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}

static void freeBUItvinfo_rdjce(GtBUItvinfo_rdjce *ptr,
                          unsigned long allocated,
                          BUstate_rdjce *state)
{
  unsigned long idx;

  for (idx=0; idx<allocated; idx++)
  {
    freeBUinfo_rdjce(&ptr[idx].info,state);
  }
  gt_free(ptr);
}

int gt_esa_bottomup_rdjce(Sequentialsuffixarrayreader *ssar,
                    BUstate_rdjce *bustate,
                    GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix = 0,
                idx,
                nonspecials,
                allocatedItvinfo = 0,
                nextfreeItvinfo = 0,
                lastsuftabvalue = 0;
  GtBUItvinfo_rdjce *lastinterval = NULL, *stackspace = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;

  PUSH_ESA_BOTTOMUP_rdjce(0,0);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  for (idx = 0; idx < nonspecials; idx++)
  {
    NEXTSEQUENTIALLCPTABVALUEWITHLAST(lcpvalue,lastsuftabvalue,ssar);
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_rdjce.lcp)
    {
      if (TOP_ESA_BOTTOMUP_rdjce.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_rdjce(firstedge,
                          TOP_ESA_BOTTOMUP_rdjce.lcp,
                          TOP_ESA_BOTTOMUP_rdjce.lb,
                          &TOP_ESA_BOTTOMUP_rdjce.info,
                          previoussuffix,bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_rdjce.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_rdjce;
      lastinterval->rb = idx;
      if (processlcpinterval_rdjce(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
        break;
      }
      if (lcpvalue <= TOP_ESA_BOTTOMUP_rdjce.lcp)
      {
        if (TOP_ESA_BOTTOMUP_rdjce.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
        if (processbranchingedge_rdjce(firstedge,
                                 TOP_ESA_BOTTOMUP_rdjce.lcp,
                                 TOP_ESA_BOTTOMUP_rdjce.lb,
                                 &TOP_ESA_BOTTOMUP_rdjce.info,
                                 lastinterval->lcp,
                                 lastinterval->lb,
                                 lastinterval->rb,
                                 &lastinterval->info,
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
    if (lcpvalue > TOP_ESA_BOTTOMUP_rdjce.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervallb = lastinterval->lb,
                      lastintervalrb = lastinterval->rb;
        PUSH_ESA_BOTTOMUP_rdjce(lcpvalue,lastintervallb);
        if (processbranchingedge_rdjce(true,
                                 TOP_ESA_BOTTOMUP_rdjce.lcp,
                                 TOP_ESA_BOTTOMUP_rdjce.lb,
                                 &TOP_ESA_BOTTOMUP_rdjce.info,
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
        PUSH_ESA_BOTTOMUP_rdjce(lcpvalue,idx);
        if (processleafedge_rdjce(true,
                            TOP_ESA_BOTTOMUP_rdjce.lcp,
                            TOP_ESA_BOTTOMUP_rdjce.lb,
                            &TOP_ESA_BOTTOMUP_rdjce.info,
                            previoussuffix,bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  gt_assert(nextfreeItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_rdjce.lcp > 0)
  {
    if (processleafedge_rdjce(false,
                        TOP_ESA_BOTTOMUP_rdjce.lcp,
                        TOP_ESA_BOTTOMUP_rdjce.lb,
                        &TOP_ESA_BOTTOMUP_rdjce.info,
                        lastsuftabvalue,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_rdjce.rb = idx;
      if (processlcpinterval_rdjce(TOP_ESA_BOTTOMUP_rdjce.lcp,
                             TOP_ESA_BOTTOMUP_rdjce.lb,
                             TOP_ESA_BOTTOMUP_rdjce.rb,
                             &TOP_ESA_BOTTOMUP_rdjce.info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }
    }
  }
  freeBUItvinfo_rdjce(stackspace, allocatedItvinfo, bustate);
  return haserr ? -1 : 0;
}
