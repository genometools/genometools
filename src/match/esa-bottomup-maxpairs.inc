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

static BUinfo_Maxpairs *allocateBUinfo_maxpairs(BUstate_Maxpairs *);
static void freeBUinfo_maxpairs(BUinfo_Maxpairs *,BUstate_Maxpairs *);
static int processleafedge_maxpairs(bool,
    unsigned long,
    unsigned long,
    BUinfo_Maxpairs *,
    unsigned long,
    BUstate_Maxpairs *,
    GtError *err);
static int processbranchingedge_maxpairs(bool firstsucc,
    unsigned long,
    unsigned long,
    BUinfo_Maxpairs *,
    unsigned long,
    unsigned long,
    unsigned long,
    BUinfo_Maxpairs *,
    BUstate_Maxpairs *,
    GtError *);

/* if processcompletenode_maxpairs is required, define the following
   macro and undefine it after the inclusion */
#ifdef USE_PROCESSCOMPLETENODE_ESA_BOTTOMUP
static int processcompletenode_maxpairs(unsigned long,
    unsigned long,
    unsigned long,
    BUinfo_Maxpairs *);
#endif

#define TOP_ESA_BOTTOMUP_MAXPAIRS   stackspace[nextfreeItvinfo-1]
#define POP_ESA_BOTTOMUP_MAXPAIRS   (stackspace + (--nextfreeItvinfo))

#define PUSH_ESA_BOTTOMUP_MAXPAIRS(LCP,LB)\
        if (nextfreeItvinfo >= allocatedItvinfo)\
        {\
          gt_assert(nextfreeItvinfo == allocatedItvinfo);\
          stackspace = allocateBUstack_maxpairs(stackspace,\
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
  BUinfo_Maxpairs *info;
} GtBUItvinfoMaxpairs;

#ifdef SKDEBUG
static void showstack_maxpairs(const GtBUItvinfoMaxpairs *stackspace,
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

static GtBUItvinfoMaxpairs *allocateBUstack_maxpairs(GtBUItvinfoMaxpairs *ptr,
                                   unsigned long currentallocated,
                                   unsigned long allocated,
                                   BUstate_Maxpairs *state)
{
  unsigned long idx;
  GtBUItvinfoMaxpairs *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof(*itvinfo) * allocated);
  gt_assert(allocated > currentallocated);
  for (idx=currentallocated; idx<allocated; idx++)
  {
    itvinfo[idx].info = allocateBUinfo_maxpairs(state);
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}

static void freeBUItvinfo_maxpairs(GtBUItvinfoMaxpairs *ptr,
                          unsigned long allocated,
                          BUstate_Maxpairs *state)
{
  unsigned long idx;

  for (idx=0; idx<allocated; idx++)
  {
    freeBUinfo_maxpairs(ptr[idx].info,state);
  }
  gt_free(ptr);
}

int gt_esa_bottomup_maxpairs(Sequentialsuffixarrayreader *ssar,
                    BUstate_Maxpairs *bustate,
                    GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix,
                idx,
                nonspecials,
                allocatedItvinfo = 0,
                nextfreeItvinfo = 0;
  GtBUItvinfoMaxpairs *lastinterval = NULL, *stackspace = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;
  int retval;

  PUSH_ESA_BOTTOMUP_MAXPAIRS(0,0);
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
    if (lcpvalue <= TOP_ESA_BOTTOMUP_MAXPAIRS.lcp)
    {
      if (TOP_ESA_BOTTOMUP_MAXPAIRS.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_maxpairs(firstedge,
                          TOP_ESA_BOTTOMUP_MAXPAIRS.lcp,
                          TOP_ESA_BOTTOMUP_MAXPAIRS.lb,
                          TOP_ESA_BOTTOMUP_MAXPAIRS.info,
                          previoussuffix,bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_MAXPAIRS.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_MAXPAIRS;
      lastinterval->rb = idx;
#ifdef USE_PROCESSCOMPLETENODE_ESA_BOTTOMUP
      processcompletenode_maxpairs(lastinterval->lcp,
                          lastinterval->lb,
                          lastinterval->rb,
                          lastinterval->info);
#endif
      if (lcpvalue <= TOP_ESA_BOTTOMUP_MAXPAIRS.lcp)
      {
        gt_assert(lastinterval->info == NULL ||
                  lastinterval->info != TOP_ESA_BOTTOMUP_MAXPAIRS.info);
        if (TOP_ESA_BOTTOMUP_MAXPAIRS.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
        if (processbranchingedge_maxpairs(firstedge,
                                 TOP_ESA_BOTTOMUP_MAXPAIRS.lcp,
                                 TOP_ESA_BOTTOMUP_MAXPAIRS.lb,
                                 TOP_ESA_BOTTOMUP_MAXPAIRS.info,
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
    if (lcpvalue > TOP_ESA_BOTTOMUP_MAXPAIRS.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervallb = lastinterval->lb,
                      lastintervalrb = lastinterval->rb;
        PUSH_ESA_BOTTOMUP_MAXPAIRS(lcpvalue,lastintervallb);
        if (processbranchingedge_maxpairs(true,
                                 TOP_ESA_BOTTOMUP_MAXPAIRS.lcp,
                                 TOP_ESA_BOTTOMUP_MAXPAIRS.lb,
                                 TOP_ESA_BOTTOMUP_MAXPAIRS.info,
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
        PUSH_ESA_BOTTOMUP_MAXPAIRS(lcpvalue,idx);
        if (processleafedge_maxpairs(true,
                            TOP_ESA_BOTTOMUP_MAXPAIRS.lcp,
                            TOP_ESA_BOTTOMUP_MAXPAIRS.lb,
                            TOP_ESA_BOTTOMUP_MAXPAIRS.info,
                            previoussuffix,bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  freeBUItvinfo_maxpairs(stackspace, allocatedItvinfo, bustate);
  return haserr ? -1 : 0;
}
