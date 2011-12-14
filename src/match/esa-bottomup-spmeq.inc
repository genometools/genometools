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

static void initBUinfo_spmeq(GtBUinfo_spmeq *,GtBUstate_spmeq *);

static void freeBUinfo_spmeq(GtBUinfo_spmeq *,GtBUstate_spmeq *);

static int processleafedge_spmeq(bool,
    unsigned long,
    unsigned long,
    GtBUinfo_spmeq *,
#ifdef GT_ESA_BOTTOM_UP_RAM
    unsigned long, /* seqnum */
    unsigned long, /* relpos */
#else
    unsigned long, /* position */
#endif
    GtBUstate_spmeq *,
    GtError *err);

#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
static int processbranchingedge_spmeq(bool firstsucc,
    unsigned long,
    unsigned long,
    GtBUinfo_spmeq *,
    unsigned long,
    unsigned long,
    unsigned long,
    GtBUinfo_spmeq *,
    GtBUstate_spmeq *,
    GtError *);
#endif

static int processlcpinterval_spmeq(unsigned long,
    unsigned long,
    unsigned long,
    GtBUinfo_spmeq *,
    GtBUstate_spmeq *,
    GtError *err);

#define TOP_ESA_BOTTOMUP_spmeq\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo-1]
#define POP_ESA_BOTTOMUP_spmeq\
        (stack->spaceGtBUItvinfo + (--stack->nextfreeGtBUItvinfo))

#define PUSH_ESA_BOTTOMUP_spmeq(LCP,LB)\
        if (stack->nextfreeGtBUItvinfo >= stack->allocatedGtBUItvinfo)\
        {\
          gt_assert(stack->nextfreeGtBUItvinfo == stack->allocatedGtBUItvinfo);\
          stack->spaceGtBUItvinfo\
            = allocateBUstack_spmeq(stack->spaceGtBUItvinfo,\
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
  GtBUinfo_spmeq info;
} GtBUItvinfo_spmeq;

typedef struct
{
  GtBUItvinfo_spmeq *spaceGtBUItvinfo;
  unsigned long allocatedGtBUItvinfo,
                nextfreeGtBUItvinfo;
} GtArrayGtBUItvinfo_spmeq;

GtArrayGtBUItvinfo_spmeq *gt_GtArrayGtBUItvinfo_new_spmeq(void)
{
  GtArrayGtBUItvinfo_spmeq *stack = gt_malloc(sizeof (*stack));

  stack->spaceGtBUItvinfo = NULL;
  stack->allocatedGtBUItvinfo = stack->nextfreeGtBUItvinfo = 0;
  return stack;
}

void gt_GtArrayGtBUItvinfo_delete_spmeq(GtArrayGtBUItvinfo_spmeq *stack,
                                  GtBUstate_spmeq *state)
{
  unsigned long idx;

  for (idx=0; idx<stack->allocatedGtBUItvinfo; idx++)
  {
    freeBUinfo_spmeq(&stack->spaceGtBUItvinfo[idx].info,state);
  }
  gt_free(stack->spaceGtBUItvinfo);
  gt_free(stack);
}

static GtBUItvinfo_spmeq *allocateBUstack_spmeq(GtBUItvinfo_spmeq *ptr,
                                   unsigned long currentallocated,
                                   unsigned long allocated,
                                   GtBUstate_spmeq *state)
{
  unsigned long idx;
  GtBUItvinfo_spmeq *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof(*itvinfo) * allocated);
  gt_assert(allocated > currentallocated);
  for (idx=currentallocated; idx<allocated; idx++)
  {
    initBUinfo_spmeq(&itvinfo[idx].info,state);
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}

#ifndef GT_ESA_BOTTOM_UP_RAM

int gt_esa_bottomup_spmeq(Sequentialsuffixarrayreader *ssar,
                    GtBUstate_spmeq *bustate,
                    GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix = 0,
                idx,
                nonspecials,
                lastsuftabvalue = 0;
  GtBUItvinfo_spmeq *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;
  GtArrayGtBUItvinfo_spmeq *stack;

  stack = gt_GtArrayGtBUItvinfo_new_spmeq();
  PUSH_ESA_BOTTOMUP_spmeq(0,0);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  for (idx = 0; idx < nonspecials; idx++)
  {
    NEXTSEQUENTIALLCPTABVALUEWITHLAST(lcpvalue,lastsuftabvalue,ssar);
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_spmeq.lcp)
    {
      if (TOP_ESA_BOTTOMUP_spmeq.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_spmeq(firstedge,
                          TOP_ESA_BOTTOMUP_spmeq.lcp,
                          TOP_ESA_BOTTOMUP_spmeq.lb,
                          &TOP_ESA_BOTTOMUP_spmeq.info,
                          previoussuffix,bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_spmeq.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_spmeq;
      lastinterval->rb = idx;
      if (processlcpinterval_spmeq != NULL &&
          processlcpinterval_spmeq(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
        break;
      }
      if (lcpvalue <= TOP_ESA_BOTTOMUP_spmeq.lcp)
      {
        if (TOP_ESA_BOTTOMUP_spmeq.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_spmeq != NULL &&
            processbranchingedge_spmeq(firstedge,
                                 TOP_ESA_BOTTOMUP_spmeq.lcp,
                                 TOP_ESA_BOTTOMUP_spmeq.lb,
                                 &TOP_ESA_BOTTOMUP_spmeq.info,
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
    if (lcpvalue > TOP_ESA_BOTTOMUP_spmeq.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallb = lastinterval->lb;
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervalrb = lastinterval->rb;
#endif
        PUSH_ESA_BOTTOMUP_spmeq(lcpvalue,lastintervallb);
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_spmeq != NULL &&
            processbranchingedge_spmeq(true,
                                 TOP_ESA_BOTTOMUP_spmeq.lcp,
                                 TOP_ESA_BOTTOMUP_spmeq.lb,
                                 &TOP_ESA_BOTTOMUP_spmeq.info,
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
        PUSH_ESA_BOTTOMUP_spmeq(lcpvalue,idx);
        if (processleafedge_spmeq(true,
                            TOP_ESA_BOTTOMUP_spmeq.lcp,
                            TOP_ESA_BOTTOMUP_spmeq.lb,
                            &TOP_ESA_BOTTOMUP_spmeq.info,
                            previoussuffix,bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_spmeq.lcp > 0)
  {
    if (processleafedge_spmeq(false,
                        TOP_ESA_BOTTOMUP_spmeq.lcp,
                        TOP_ESA_BOTTOMUP_spmeq.lb,
                        &TOP_ESA_BOTTOMUP_spmeq.info,
                        lastsuftabvalue,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_spmeq.rb = idx;
      if (processlcpinterval_spmeq(TOP_ESA_BOTTOMUP_spmeq.lcp,
                             TOP_ESA_BOTTOMUP_spmeq.lb,
                             TOP_ESA_BOTTOMUP_spmeq.rb,
                             &TOP_ESA_BOTTOMUP_spmeq.info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_GtArrayGtBUItvinfo_delete_spmeq(stack,bustate);
  return haserr ? -1 : 0;
}

#else

int gt_esa_bottomup_RAM_spmeq(const unsigned long *seqnum_relpos_bucket,
                        const GtSeqnumrelpos *snrp,
                        const uint16_t *lcptab_bucket,
                        unsigned long nonspecials,
                        GtArrayGtBUItvinfo_spmeq *stack,
                        GtBUstate_spmeq *bustate,
                        GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                tmp,
                previousseqnum = 0,
                previousrelpos = 0,
                idx;
  GtBUItvinfo_spmeq *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;

  gt_assert(nonspecials > 0);
  PUSH_ESA_BOTTOMUP_spmeq(0,0);
  for (idx = 0; idx < nonspecials-1; idx++)
  {
    lcpvalue = (unsigned long) lcptab_bucket[idx+1];
    tmp = seqnum_relpos_bucket[idx];
    previousseqnum = gt_seqnumrelpos_decode_seqnum(snrp,tmp);
    previousrelpos = gt_seqnumrelpos_decode_relpos(snrp,tmp);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_spmeq.lcp)
    {
      if (TOP_ESA_BOTTOMUP_spmeq.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_spmeq(firstedge,
                          TOP_ESA_BOTTOMUP_spmeq.lcp,
                          TOP_ESA_BOTTOMUP_spmeq.lb,
                          &TOP_ESA_BOTTOMUP_spmeq.info,
                          previousseqnum,previousrelpos,bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_spmeq.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_spmeq;
      lastinterval->rb = idx;
      if (processlcpinterval_spmeq != NULL &&
          processlcpinterval_spmeq(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
        break;
      }
      if (lcpvalue <= TOP_ESA_BOTTOMUP_spmeq.lcp)
      {
        if (TOP_ESA_BOTTOMUP_spmeq.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_spmeq != NULL &&
            processbranchingedge_spmeq(firstedge,
                                 TOP_ESA_BOTTOMUP_spmeq.lcp,
                                 TOP_ESA_BOTTOMUP_spmeq.lb,
                                 &TOP_ESA_BOTTOMUP_spmeq.info,
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
    if (lcpvalue > TOP_ESA_BOTTOMUP_spmeq.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallb = lastinterval->lb;
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervalrb = lastinterval->rb;
#endif
        PUSH_ESA_BOTTOMUP_spmeq(lcpvalue,lastintervallb);
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_spmeq != NULL &&
            processbranchingedge_spmeq(true,
                                 TOP_ESA_BOTTOMUP_spmeq.lcp,
                                 TOP_ESA_BOTTOMUP_spmeq.lb,
                                 &TOP_ESA_BOTTOMUP_spmeq.info,
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
        PUSH_ESA_BOTTOMUP_spmeq(lcpvalue,idx);
        if (processleafedge_spmeq(true,
                            TOP_ESA_BOTTOMUP_spmeq.lcp,
                            TOP_ESA_BOTTOMUP_spmeq.lb,
                            &TOP_ESA_BOTTOMUP_spmeq.info,
                            previousseqnum,previousrelpos,bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_spmeq.lcp > 0)
  {
    unsigned long lastseqnum, lastrelpos;

    tmp = seqnum_relpos_bucket[nonspecials-1];
    lastseqnum = gt_seqnumrelpos_decode_seqnum(snrp,tmp);
    lastrelpos = gt_seqnumrelpos_decode_relpos(snrp,tmp);
    if (processleafedge_spmeq(false,
                        TOP_ESA_BOTTOMUP_spmeq.lcp,
                        TOP_ESA_BOTTOMUP_spmeq.lb,
                        &TOP_ESA_BOTTOMUP_spmeq.info,
                        lastseqnum,lastrelpos,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_spmeq.rb = idx;
      if (processlcpinterval_spmeq(TOP_ESA_BOTTOMUP_spmeq.lcp,
                             TOP_ESA_BOTTOMUP_spmeq.lb,
                             TOP_ESA_BOTTOMUP_spmeq.rb,
                             &TOP_ESA_BOTTOMUP_spmeq.info,
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
