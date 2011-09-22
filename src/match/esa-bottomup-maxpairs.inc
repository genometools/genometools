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
#ifdef GT_ESA_BOTTOM_UP_RAM
#include "seqnumrelpos.h"
#endif

static void initBUinfo_maxpairs(GtBUinfo_maxpairs *,GtBUstate_maxpairs *);

static void freeBUinfo_maxpairs(GtBUinfo_maxpairs *,GtBUstate_maxpairs *);

static int processleafedge_maxpairs(bool,
    unsigned long,
    unsigned long,
    GtBUinfo_maxpairs *,
    unsigned long,
#ifdef GT_ESA_BOTTOM_UP_RAM
    unsigned long,
    unsigned long,
#endif
    GtBUstate_maxpairs *,
    GtError *err);

#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
static int processbranchingedge_maxpairs(bool firstsucc,
    unsigned long,
    unsigned long,
    GtBUinfo_maxpairs *,
    unsigned long,
    unsigned long,
    unsigned long,
    GtBUinfo_maxpairs *,
    GtBUstate_maxpairs *,
    GtError *);
#endif

static int processlcpinterval_maxpairs(unsigned long,
    unsigned long,
    unsigned long,
    GtBUinfo_maxpairs *,
    GtBUstate_maxpairs *,
    GtError *err);

#define TOP_ESA_BOTTOMUP_maxpairs\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo-1]
#define POP_ESA_BOTTOMUP_maxpairs\
        (stack->spaceGtBUItvinfo + (--stack->nextfreeGtBUItvinfo))

#define PUSH_ESA_BOTTOMUP_maxpairs(LCP,LB)\
        if (stack->nextfreeGtBUItvinfo >= stack->allocatedGtBUItvinfo)\
        {\
          gt_assert(stack->nextfreeGtBUItvinfo == stack->allocatedGtBUItvinfo);\
          stack->spaceGtBUItvinfo\
            = allocateBUstack_maxpairs(stack->spaceGtBUItvinfo,\
                              stack->allocatedGtBUItvinfo,\
                              stack->allocatedGtBUItvinfo+incrementstacksize,\
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
  GtBUinfo_maxpairs info;
} GtBUItvinfo_maxpairs;

typedef struct
{
  GtBUItvinfo_maxpairs *spaceGtBUItvinfo;
  unsigned long allocatedGtBUItvinfo,
                nextfreeGtBUItvinfo;
} GtArrayGtBUItvinfo_maxpairs;

GtArrayGtBUItvinfo_maxpairs *gt_GtArrayGtBUItvinfo_new_maxpairs(void)
{
  GtArrayGtBUItvinfo_maxpairs *stack = gt_malloc(sizeof (*stack));

  stack->spaceGtBUItvinfo = NULL;
  stack->allocatedGtBUItvinfo = stack->nextfreeGtBUItvinfo = 0;
  return stack;
}

void gt_GtArrayGtBUItvinfo_delete_maxpairs(GtArrayGtBUItvinfo_maxpairs *stack,
                                  GtBUstate_maxpairs *state)
{
  unsigned long idx;

  for (idx=0; idx<stack->allocatedGtBUItvinfo; idx++)
  {
    freeBUinfo_maxpairs(&stack->spaceGtBUItvinfo[idx].info,state);
  }
  gt_free(stack->spaceGtBUItvinfo);
  gt_free(stack);
}

static GtBUItvinfo_maxpairs *allocateBUstack_maxpairs(GtBUItvinfo_maxpairs *ptr,
                                   unsigned long currentallocated,
                                   unsigned long allocated,
                                   GtBUstate_maxpairs *state)
{
  unsigned long idx;
  GtBUItvinfo_maxpairs *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof(*itvinfo) * allocated);
  gt_assert(allocated > currentallocated);
  for (idx=currentallocated; idx<allocated; idx++)
  {
    initBUinfo_maxpairs(&itvinfo[idx].info,state);
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}

#ifndef GT_ESA_BOTTOM_UP_RAM

int gt_esa_bottomup_maxpairs(Sequentialsuffixarrayreader *ssar,
                    GtBUstate_maxpairs *bustate,
                    GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix = 0,
                idx,
                nonspecials,
                lastsuftabvalue = 0;
  GtBUItvinfo_maxpairs *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;
  GtArrayGtBUItvinfo_maxpairs *stack;

  stack = gt_GtArrayGtBUItvinfo_new_maxpairs();
  PUSH_ESA_BOTTOMUP_maxpairs(0,0);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  for (idx = 0; idx < nonspecials; idx++)
  {
    NEXTSEQUENTIALLCPTABVALUEWITHLAST(lcpvalue,lastsuftabvalue,ssar);
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_maxpairs.lcp)
    {
      if (TOP_ESA_BOTTOMUP_maxpairs.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_maxpairs(firstedge,
                          TOP_ESA_BOTTOMUP_maxpairs.lcp,
                          TOP_ESA_BOTTOMUP_maxpairs.lb,
                          &TOP_ESA_BOTTOMUP_maxpairs.info,
                          previoussuffix,bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_maxpairs.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_maxpairs;
      lastinterval->rb = idx;
      if (processlcpinterval_maxpairs != NULL &&
          processlcpinterval_maxpairs(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
        break;
      }
      if (lcpvalue <= TOP_ESA_BOTTOMUP_maxpairs.lcp)
      {
        if (TOP_ESA_BOTTOMUP_maxpairs.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_maxpairs != NULL &&
            processbranchingedge_maxpairs(firstedge,
                                 TOP_ESA_BOTTOMUP_maxpairs.lcp,
                                 TOP_ESA_BOTTOMUP_maxpairs.lb,
                                 &TOP_ESA_BOTTOMUP_maxpairs.info,
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
#endif
        lastinterval = NULL;
      }
    }
    if (haserr)
    {
      break;
    }
    if (lcpvalue > TOP_ESA_BOTTOMUP_maxpairs.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallb = lastinterval->lb;
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervalrb = lastinterval->rb;
#endif
        PUSH_ESA_BOTTOMUP_maxpairs(lcpvalue,lastintervallb);
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_maxpairs != NULL &&
            processbranchingedge_maxpairs(true,
                                 TOP_ESA_BOTTOMUP_maxpairs.lcp,
                                 TOP_ESA_BOTTOMUP_maxpairs.lb,
                                 &TOP_ESA_BOTTOMUP_maxpairs.info,
                                 lastintervallcp,
                                 lastintervallb,
                                 lastintervalrb,
                                 NULL,
                                 bustate,err) != 0)
        {
          haserr = true;
          break;
        }
#endif
        lastinterval = NULL;
      } else
      {
        PUSH_ESA_BOTTOMUP_maxpairs(lcpvalue,idx);
        if (processleafedge_maxpairs(true,
                            TOP_ESA_BOTTOMUP_maxpairs.lcp,
                            TOP_ESA_BOTTOMUP_maxpairs.lb,
                            &TOP_ESA_BOTTOMUP_maxpairs.info,
                            previoussuffix,bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_maxpairs.lcp > 0)
  {
    if (processleafedge_maxpairs(false,
                        TOP_ESA_BOTTOMUP_maxpairs.lcp,
                        TOP_ESA_BOTTOMUP_maxpairs.lb,
                        &TOP_ESA_BOTTOMUP_maxpairs.info,
                        lastsuftabvalue,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_maxpairs.rb = idx;
      if (processlcpinterval_maxpairs(TOP_ESA_BOTTOMUP_maxpairs.lcp,
                             TOP_ESA_BOTTOMUP_maxpairs.lb,
                             TOP_ESA_BOTTOMUP_maxpairs.rb,
                             &TOP_ESA_BOTTOMUP_maxpairs.info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_GtArrayGtBUItvinfo_delete_maxpairs(stack,bustate);
  return haserr ? -1 : 0;
}

#else

int gt_esa_bottomup_RAM_maxpairs(const unsigned long *suftab_bucket,
                        GT_UNUSED const GtSeqnumrelpostab *snrp,
                        const uint16_t *lcptab_bucket,
                        unsigned long nonspecials,
                        GtArrayGtBUItvinfo_maxpairs *stack,
                        GtBUstate_maxpairs *bustate,
                        GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix = 0,
                previousseqnum = 0,
                previousrelpos = 0,
                idx;
  GtBUItvinfo_maxpairs *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;

  gt_assert(nonspecials > 0);
  PUSH_ESA_BOTTOMUP_maxpairs(0,0);
  for (idx = 0; idx < nonspecials-1; idx++)
  {
    lcpvalue = (unsigned long) lcptab_bucket[idx+1];
    previoussuffix = suftab_bucket[idx];
    previousseqnum = gt_seqnumrelpostab_seqnum(snrp,idx);
    previousrelpos = gt_seqnumrelpostab_relpos(snrp,idx);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_maxpairs.lcp)
    {
      if (TOP_ESA_BOTTOMUP_maxpairs.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_maxpairs(firstedge,
                          TOP_ESA_BOTTOMUP_maxpairs.lcp,
                          TOP_ESA_BOTTOMUP_maxpairs.lb,
                          &TOP_ESA_BOTTOMUP_maxpairs.info,
                          previoussuffix,previousseqnum,previousrelpos,
                          bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_maxpairs.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_maxpairs;
      lastinterval->rb = idx;
      if (processlcpinterval_maxpairs != NULL &&
          processlcpinterval_maxpairs(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
        break;
      }
      if (lcpvalue <= TOP_ESA_BOTTOMUP_maxpairs.lcp)
      {
        if (TOP_ESA_BOTTOMUP_maxpairs.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_maxpairs != NULL &&
            processbranchingedge_maxpairs(firstedge,
                                 TOP_ESA_BOTTOMUP_maxpairs.lcp,
                                 TOP_ESA_BOTTOMUP_maxpairs.lb,
                                 &TOP_ESA_BOTTOMUP_maxpairs.info,
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
#endif
        lastinterval = NULL;
      }
    }
    if (haserr)
    {
      break;
    }
    if (lcpvalue > TOP_ESA_BOTTOMUP_maxpairs.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallb = lastinterval->lb;
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervalrb = lastinterval->rb;
#endif
        PUSH_ESA_BOTTOMUP_maxpairs(lcpvalue,lastintervallb);
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_maxpairs != NULL &&
            processbranchingedge_maxpairs(true,
                                 TOP_ESA_BOTTOMUP_maxpairs.lcp,
                                 TOP_ESA_BOTTOMUP_maxpairs.lb,
                                 &TOP_ESA_BOTTOMUP_maxpairs.info,
                                 lastintervallcp,
                                 lastintervallb,
                                 lastintervalrb,
                                 NULL,
                                 bustate,err) != 0)
        {
          haserr = true;
          break;
        }
#endif
        lastinterval = NULL;
      } else
      {
        PUSH_ESA_BOTTOMUP_maxpairs(lcpvalue,idx);
        if (processleafedge_maxpairs(true,
                            TOP_ESA_BOTTOMUP_maxpairs.lcp,
                            TOP_ESA_BOTTOMUP_maxpairs.lb,
                            &TOP_ESA_BOTTOMUP_maxpairs.info,
                            previoussuffix,previousseqnum,previousrelpos,
                            bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_maxpairs.lcp > 0)
  {
    unsigned long lastsuftabvalue, lastseqnum, lastrelpos;

    lastsuftabvalue = suftab_bucket[nonspecials-1];
    lastseqnum = gt_seqnumrelpostab_seqnum(snrp,nonspecials-1);
    lastrelpos = gt_seqnumrelpostab_relpos(snrp,nonspecials-1);
    if (processleafedge_maxpairs(false,
                        TOP_ESA_BOTTOMUP_maxpairs.lcp,
                        TOP_ESA_BOTTOMUP_maxpairs.lb,
                        &TOP_ESA_BOTTOMUP_maxpairs.info,
                        lastsuftabvalue,lastseqnum,lastrelpos,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_maxpairs.rb = idx;
      if (processlcpinterval_maxpairs(TOP_ESA_BOTTOMUP_maxpairs.lcp,
                             TOP_ESA_BOTTOMUP_maxpairs.lb,
                             TOP_ESA_BOTTOMUP_maxpairs.rb,
                             &TOP_ESA_BOTTOMUP_maxpairs.info,
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

#endif
