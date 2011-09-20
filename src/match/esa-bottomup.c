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

#define TOP_ESA_BOTTOMUP\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo-1]
#define POP_ESA_BOTTOMUP\
        (stack->spaceGtBUItvinfo + (--stack->nextfreeGtBUItvinfo))

#define PUSH_ESA_BOTTOMUP(LCP,LB)\
        if (stack->nextfreeGtBUItvinfo >= stack->allocatedGtBUItvinfo)\
        {\
          gt_assert(stack->nextfreeGtBUItvinfo == stack->allocatedGtBUItvinfo);\
          stack->spaceGtBUItvinfo\
            = allocateBUstack(stack->spaceGtBUItvinfo,\
                              stack->allocatedGtBUItvinfo,\
                              stack->allocatedGtBUItvinfo+incrementstacksize,\
                              allocateBUinfo,\
                              bustate);\
          stack->allocatedGtBUItvinfo += incrementstacksize;\
        }\
        gt_assert(stack->spaceGtBUItvinfo != NULL);\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo].lcp = LCP;\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo].lb = LB;\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo].rb = ULONG_MAX;\
        stack->nextfreeGtBUItvinfo++

typedef struct
{
  unsigned long lcp, lb, rb;
  GtBUinfo *info;
} GtBUItvinfo;

struct GtArrayGtBUItvinfo
{
  GtBUItvinfo *spaceGtBUItvinfo;
  unsigned long allocatedGtBUItvinfo,
                nextfreeGtBUItvinfo;
};

GtArrayGtBUItvinfo *gt_GtArrayGtBUItvinfo_new(void)
{
  GtArrayGtBUItvinfo *stack = gt_malloc(sizeof (*stack));

  stack->spaceGtBUItvinfo = NULL;
  stack->allocatedGtBUItvinfo = stack->nextfreeGtBUItvinfo = 0;
  return stack;
}

void gt_GtArrayGtBUItvinfo_delete(GtArrayGtBUItvinfo *stack,
                                  void (*freeBUinfo)(GtBUinfo *,GtBUstate *),
                                  GtBUstate *state)
{
  unsigned long idx;

  for (idx=0; idx<stack->allocatedGtBUItvinfo; idx++)
  {
    if (freeBUinfo != NULL)
    {
      freeBUinfo(stack->spaceGtBUItvinfo[idx].info,state);
    }
  }
  gt_free(stack->spaceGtBUItvinfo);
  gt_free(stack);
}

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
                    int (*processlcpinterval)(unsigned long,
                                              unsigned long,
                                              unsigned long,
                                              GtBUinfo *,
                                              GtBUstate *,
                                              GtError *err),
                    GtBUstate *bustate,
                    GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix = 0,
                idx,
                nonspecials,
                lastsuftabvalue = 0;
  GtBUItvinfo *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;
  GtArrayGtBUItvinfo *stack;

  stack = gt_GtArrayGtBUItvinfo_new();
  PUSH_ESA_BOTTOMUP(0,0);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  for (idx = 0; idx < nonspecials; idx++)
  {
    NEXTSEQUENTIALLCPTABVALUEWITHLAST(lcpvalue,lastsuftabvalue,ssar);
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
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
      if (processlcpinterval != NULL &&
          processlcpinterval(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             lastinterval->info,
                             bustate,
                             err) != 0)
      {
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
        if (processbranchingedge != NULL &&
            processbranchingedge(firstedge,
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
        if (processbranchingedge != NULL &&
            processbranchingedge(true,
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
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP.lcp > 0)
  {
    if (processleafedge(false,
                        TOP_ESA_BOTTOMUP.lcp,
                        TOP_ESA_BOTTOMUP.lb,
                        TOP_ESA_BOTTOMUP.info,
                        lastsuftabvalue,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP.rb = idx;
      if (processlcpinterval(TOP_ESA_BOTTOMUP.lcp,
                             TOP_ESA_BOTTOMUP.lb,
                             TOP_ESA_BOTTOMUP.rb,
                             TOP_ESA_BOTTOMUP.info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_GtArrayGtBUItvinfo_delete(stack,freeBUinfo,bustate);
  return haserr ? -1 : 0;
}

int gt_esa_bottomup_RAM(const unsigned long *suftab,
                        const uint16_t *lcptab_bucket,
                        unsigned long nonspecials,
                        GtArrayGtBUItvinfo *stack,
                        GtBUinfo *(*allocateBUinfo)(GtBUstate *),
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
                        int (*processlcpinterval)(unsigned long,
                                                  unsigned long,
                                                  unsigned long,
                                                  GtBUinfo *,
                                                  GtBUstate *,
                                                  GtError *err),
                        GtBUstate *bustate,
                        GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix = 0,
                idx;
  GtBUItvinfo *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;

  gt_assert(nonspecials > 0);
  PUSH_ESA_BOTTOMUP(0,0);
  for (idx = 0; idx < nonspecials-1; idx++)
  {
    lcpvalue = (unsigned long) lcptab_bucket[idx+1];
    previoussuffix = suftab[idx];
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
      if (processlcpinterval != NULL &&
          processlcpinterval(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             lastinterval->info,
                             bustate,
                             err) != 0)
      {
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
        if (processbranchingedge != NULL &&
            processbranchingedge(firstedge,
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
        if (processbranchingedge != NULL &&
            processbranchingedge(true,
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
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP.lcp > 0)
  {
    unsigned long lastsuftabvalue = suftab[nonspecials-1];
    if (processleafedge(false,
                        TOP_ESA_BOTTOMUP.lcp,
                        TOP_ESA_BOTTOMUP.lb,
                        TOP_ESA_BOTTOMUP.info,
                        lastsuftabvalue,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP.rb = idx;
      if (processlcpinterval(TOP_ESA_BOTTOMUP.lcp,
                             TOP_ESA_BOTTOMUP.lb,
                             TOP_ESA_BOTTOMUP.rb,
                             TOP_ESA_BOTTOMUP.info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }
    }
  }
  stack->nextfreeGtBUItvinfo = 0; /* empty the stack */
  return haserr ? -1 : 0;
}
