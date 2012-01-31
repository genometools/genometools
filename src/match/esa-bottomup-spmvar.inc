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

static void initBUinfo_spmvar(GtBUinfo_spmvar *,GtBUstate_spmvar *);

static void freeBUinfo_spmvar(GtBUinfo_spmvar *,GtBUstate_spmvar *);

static int processleafedge_spmvar(bool,
    unsigned long,
    unsigned long,
    GtBUinfo_spmvar *,
#ifdef GT_ESA_BOTTOM_UP_RAM
    unsigned long, /* seqnum */
    unsigned long, /* relpos */
#else
    unsigned long, /* position */
#endif
    GtBUstate_spmvar *,
    GtError *err);

#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
static int processbranchingedge_spmvar(bool firstsucc,
    unsigned long,
    unsigned long,
    GtBUinfo_spmvar *,
    unsigned long,
    unsigned long,
    unsigned long,
    GtBUinfo_spmvar *,
    GtBUstate_spmvar *,
    GtError *);
#endif

static int processlcpinterval_spmvar(unsigned long,
    unsigned long,
    unsigned long,
    GtBUinfo_spmvar *,
    GtBUstate_spmvar *,
    GtError *err);

#define TOP_ESA_BOTTOMUP_spmvar\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo-1]
#define POP_ESA_BOTTOMUP_spmvar\
        (stack->spaceGtBUItvinfo + (--stack->nextfreeGtBUItvinfo))

#define PUSH_ESA_BOTTOMUP_spmvar(LCP,LB)\
        if (stack->nextfreeGtBUItvinfo >= stack->allocatedGtBUItvinfo)\
        {\
          gt_assert(stack->nextfreeGtBUItvinfo == stack->allocatedGtBUItvinfo);\
          stack->spaceGtBUItvinfo\
            = allocateBUstack_spmvar(stack->spaceGtBUItvinfo,\
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
  GtBUinfo_spmvar info;
} GtBUItvinfo_spmvar;

typedef struct
{
  GtBUItvinfo_spmvar *spaceGtBUItvinfo;
  unsigned long allocatedGtBUItvinfo,
                nextfreeGtBUItvinfo;
} GtArrayGtBUItvinfo_spmvar;

GtArrayGtBUItvinfo_spmvar *gt_GtArrayGtBUItvinfo_new_spmvar(void)
{
  GtArrayGtBUItvinfo_spmvar *stack = gt_malloc(sizeof (*stack));

  stack->spaceGtBUItvinfo = NULL;
  stack->allocatedGtBUItvinfo = stack->nextfreeGtBUItvinfo = 0;
  return stack;
}

void gt_GtArrayGtBUItvinfo_delete_spmvar(GtArrayGtBUItvinfo_spmvar *stack,
                                  GtBUstate_spmvar *state)
{
  unsigned long idx;

  for (idx=0; idx<stack->allocatedGtBUItvinfo; idx++)
  {
    freeBUinfo_spmvar(&stack->spaceGtBUItvinfo[idx].info,state);
  }
  gt_free(stack->spaceGtBUItvinfo);
  gt_free(stack);
}

static GtBUItvinfo_spmvar *allocateBUstack_spmvar(GtBUItvinfo_spmvar *ptr,
                                   unsigned long currentallocated,
                                   unsigned long allocated,
                                   GtBUstate_spmvar *state)
{
  unsigned long idx;
  GtBUItvinfo_spmvar *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof (*itvinfo) * allocated);
  gt_assert(allocated > currentallocated);
  for (idx=currentallocated; idx<allocated; idx++)
  {
    initBUinfo_spmvar(&itvinfo[idx].info,state);
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}

#ifndef GT_ESA_BOTTOM_UP_RAM

int gt_esa_bottomup_spmvar(Sequentialsuffixarrayreader *ssar,
                    GtBUstate_spmvar *bustate,
                    GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix = 0,
                idx,
                nonspecials,
                lastsuftabvalue = 0;
  GtBUItvinfo_spmvar *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;
  GtArrayGtBUItvinfo_spmvar *stack;

  stack = gt_GtArrayGtBUItvinfo_new_spmvar();
  PUSH_ESA_BOTTOMUP_spmvar(0,0);
  nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  for (idx = 0; idx < nonspecials; idx++)
  {
    NEXTSEQUENTIALLCPTABVALUEWITHLAST(lcpvalue,lastsuftabvalue,ssar);
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      if (TOP_ESA_BOTTOMUP_spmvar.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_spmvar(firstedge,
                          TOP_ESA_BOTTOMUP_spmvar.lcp,
                          TOP_ESA_BOTTOMUP_spmvar.lb,
                          &TOP_ESA_BOTTOMUP_spmvar.info,
                          previoussuffix,bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_spmvar;
      lastinterval->rb = idx;
      if (processlcpinterval_spmvar(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
        break;
      }
      if (lcpvalue <= TOP_ESA_BOTTOMUP_spmvar.lcp)
      {
        if (TOP_ESA_BOTTOMUP_spmvar.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_spmvar(firstedge,
                                 TOP_ESA_BOTTOMUP_spmvar.lcp,
                                 TOP_ESA_BOTTOMUP_spmvar.lb,
                                 &TOP_ESA_BOTTOMUP_spmvar.info,
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
    if (lcpvalue > TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallb = lastinterval->lb;
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervalrb = lastinterval->rb;
#endif
        PUSH_ESA_BOTTOMUP_spmvar(lcpvalue,lastintervallb);
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_spmvar(true,
                                 TOP_ESA_BOTTOMUP_spmvar.lcp,
                                 TOP_ESA_BOTTOMUP_spmvar.lb,
                                 &TOP_ESA_BOTTOMUP_spmvar.info,
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
        PUSH_ESA_BOTTOMUP_spmvar(lcpvalue,idx);
        if (processleafedge_spmvar(true,
                            TOP_ESA_BOTTOMUP_spmvar.lcp,
                            TOP_ESA_BOTTOMUP_spmvar.lb,
                            &TOP_ESA_BOTTOMUP_spmvar.info,
                            previoussuffix,bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_spmvar.lcp > 0)
  {
    if (processleafedge_spmvar(false,
                        TOP_ESA_BOTTOMUP_spmvar.lcp,
                        TOP_ESA_BOTTOMUP_spmvar.lb,
                        &TOP_ESA_BOTTOMUP_spmvar.info,
                        lastsuftabvalue,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_spmvar.rb = idx;
      if (processlcpinterval_spmvar(TOP_ESA_BOTTOMUP_spmvar.lcp,
                             TOP_ESA_BOTTOMUP_spmvar.lb,
                             TOP_ESA_BOTTOMUP_spmvar.rb,
                             &TOP_ESA_BOTTOMUP_spmvar.info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_GtArrayGtBUItvinfo_delete_spmvar(stack,bustate);
  return haserr ? -1 : 0;
}

#else

int gt_esa_bottomup_RAM_spmvar(const unsigned long *seqnum_relpos_bucket,
                        const GtSeqnumrelpos *snrp,
                        const uint16_t *lcptab_bucket,
                        unsigned long nonspecials,
                        GtArrayGtBUItvinfo_spmvar *stack,
                        GtBUstate_spmvar *bustate,
                        GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                tmp,
                previousseqnum = 0,
                previousrelpos = 0,
                idx;
  GtBUItvinfo_spmvar *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;

  gt_assert(nonspecials > 0);
  PUSH_ESA_BOTTOMUP_spmvar(0,0);
  for (idx = 0; idx < nonspecials-1; idx++)
  {
    lcpvalue = (unsigned long) lcptab_bucket[idx+1];
    tmp = seqnum_relpos_bucket[idx];
    previousseqnum = gt_seqnumrelpos_decode_seqnum(snrp,tmp);
    previousrelpos = gt_seqnumrelpos_decode_relpos(snrp,tmp);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      if (TOP_ESA_BOTTOMUP_spmvar.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_spmvar(firstedge,
                          TOP_ESA_BOTTOMUP_spmvar.lcp,
                          TOP_ESA_BOTTOMUP_spmvar.lb,
                          &TOP_ESA_BOTTOMUP_spmvar.info,
                          previousseqnum,previousrelpos,bustate,err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_spmvar;
      lastinterval->rb = idx;
      if (processlcpinterval_spmvar(lastinterval->lcp,
                             lastinterval->lb,
                             lastinterval->rb,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
        break;
      }
      if (lcpvalue <= TOP_ESA_BOTTOMUP_spmvar.lcp)
      {
        if (TOP_ESA_BOTTOMUP_spmvar.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_spmvar(firstedge,
                                 TOP_ESA_BOTTOMUP_spmvar.lcp,
                                 TOP_ESA_BOTTOMUP_spmvar.lb,
                                 &TOP_ESA_BOTTOMUP_spmvar.info,
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
    if (lcpvalue > TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallb = lastinterval->lb;
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        unsigned long lastintervallcp = lastinterval->lcp,
                      lastintervalrb = lastinterval->rb;
#endif
        PUSH_ESA_BOTTOMUP_spmvar(lcpvalue,lastintervallb);
#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
        if (processbranchingedge_spmvar(true,
                                 TOP_ESA_BOTTOMUP_spmvar.lcp,
                                 TOP_ESA_BOTTOMUP_spmvar.lb,
                                 &TOP_ESA_BOTTOMUP_spmvar.info,
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
        PUSH_ESA_BOTTOMUP_spmvar(lcpvalue,idx);
        if (processleafedge_spmvar(true,
                            TOP_ESA_BOTTOMUP_spmvar.lcp,
                            TOP_ESA_BOTTOMUP_spmvar.lb,
                            &TOP_ESA_BOTTOMUP_spmvar.info,
                            previousseqnum,previousrelpos,bustate,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_spmvar.lcp > 0)
  {
    unsigned long lastseqnum, lastrelpos;

    tmp = seqnum_relpos_bucket[nonspecials-1];
    lastseqnum = gt_seqnumrelpos_decode_seqnum(snrp,tmp);
    lastrelpos = gt_seqnumrelpos_decode_relpos(snrp,tmp);
    if (processleafedge_spmvar(false,
                        TOP_ESA_BOTTOMUP_spmvar.lcp,
                        TOP_ESA_BOTTOMUP_spmvar.lb,
                        &TOP_ESA_BOTTOMUP_spmvar.info,
                        lastseqnum,lastrelpos,bustate,err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_spmvar.rb = idx;
      if (processlcpinterval_spmvar(TOP_ESA_BOTTOMUP_spmvar.lcp,
                             TOP_ESA_BOTTOMUP_spmvar.lb,
                             TOP_ESA_BOTTOMUP_spmvar.rb,
                             &TOP_ESA_BOTTOMUP_spmvar.info,
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
/* vim: ft=c
*/
