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
    GtBUinfo_maxpairs *,
#ifdef GT_ESA_BOTTOM_UP_RAM
    unsigned long, /* seqnum */
    unsigned long, /* relpos */
#else
    unsigned long, /* position */
#endif
    GtBUstate_maxpairs *,
    GtError *err);

#ifndef GT_ESA_BOTTOM_UP_IGNORE_PROCESSBRANCHING_EDGE
static int processbranchingedge_maxpairs(bool firstsucc,
    unsigned long,
    GtBUinfo_maxpairs *,
    unsigned long,
    unsigned long,
    GtBUinfo_maxpairs *,
    GtBUstate_maxpairs *,
    GtError *);
#endif

static int processlcpinterval_maxpairs(unsigned long,
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
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo++].rb = ULONG_MAX

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

  itvinfo = gt_realloc(ptr,sizeof (*itvinfo) * allocated);
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
                numberofsuffixes,
                lastsuftabvalue = 0;
  GtBUItvinfo_maxpairs *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;
  GtArrayGtBUItvinfo_maxpairs *stack;

  stack = gt_GtArrayGtBUItvinfo_new_maxpairs();
  PUSH_ESA_BOTTOMUP_maxpairs(0,0);
  numberofsuffixes = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  for (idx = 0; idx < numberofsuffixes; idx++)
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
                          &TOP_ESA_BOTTOMUP_maxpairs.info,
                          previoussuffix,
                          bustate,
                          err) != 0)
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
      if (processlcpinterval_maxpairs(lastinterval->lcp,
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
        if (processbranchingedge_maxpairs(firstedge,
                                 TOP_ESA_BOTTOMUP_maxpairs.lcp,
                                 &TOP_ESA_BOTTOMUP_maxpairs.info,
                                 lastinterval->lcp,
                                 lastinterval->rb - lastinterval->lb + 1,
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
        if (processbranchingedge_maxpairs(true,
                                 TOP_ESA_BOTTOMUP_maxpairs.lcp,
                                 &TOP_ESA_BOTTOMUP_maxpairs.info,
                                 lastintervallcp,
                                 lastintervalrb - lastintervallb + 1,
                                 NULL,
                                 bustate,
                                 err) != 0)
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
                            &TOP_ESA_BOTTOMUP_maxpairs.info,
                            previoussuffix,
                            bustate,
                            err) != 0)
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
                        &TOP_ESA_BOTTOMUP_maxpairs.info,
                        lastsuftabvalue,
                        bustate,
                        err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_maxpairs.rb = idx;
      if (processlcpinterval_maxpairs(TOP_ESA_BOTTOMUP_maxpairs.lcp,
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

int gt_esa_bottomup_RAM_maxpairs(const unsigned long *bucketofsuffixes,
                        const GtBUlcptype_maxpairs *lcptab_bucket,
                        unsigned long numberofsuffixes,
                        GtArrayGtBUItvinfo_maxpairs *stack,
                        GtBUstate_maxpairs *bustate,
#ifdef GT_ESA_BOTTOM_UP_SEQNUM_RELPOS
                        const GtSeqnumrelpos *snrp,
#endif
                        GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix,
                lastsuftabvalue,
#ifdef GT_ESA_BOTTOM_UP_SEQNUM_RELPOS
                previousseqnum = 0,
                previousrelpos = 0,
#endif
                idx;
  GtBUItvinfo_maxpairs *lastinterval = NULL;
  bool haserr = false, 
       firstedge, 
       firstedgefromroot = true; /* must be part of state */

  gt_assert(numberofsuffixes > 0);
  PUSH_ESA_BOTTOMUP_maxpairs(0,0);
  for (idx = 0; idx < numberofsuffixes-1; idx++)
  {
    lcpvalue = (unsigned long) lcptab_bucket[idx+1];
    previoussuffix = bucketofsuffixes[idx];
#ifdef GT_ESA_BOTTOM_UP_SEQNUM_RELPOS
    previousseqnum = gt_seqnumrelpos_decode_seqnum(snrp,previoussuffix);
    previousrelpos = gt_seqnumrelpos_decode_relpos(snrp,previoussuffix);
#endif
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
                          &TOP_ESA_BOTTOMUP_maxpairs.info,
#ifdef GT_ESA_BOTTOM_UP_SEQNUM_RELPOS
                          previousseqnum,
                          previousrelpos,
#else
                          previoussuffix,
#endif
                          bustate,
                          err) != 0)
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
      if (processlcpinterval_maxpairs(lastinterval->lcp,
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
        if (processbranchingedge_maxpairs(firstedge,
                                 TOP_ESA_BOTTOMUP_maxpairs.lcp,
                                 &TOP_ESA_BOTTOMUP_maxpairs.info,
                                 lastinterval->lcp,
                                 lastinterval->rb - lastinterval->lb + 1,
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
        if (processbranchingedge_maxpairs(true,
                                 TOP_ESA_BOTTOMUP_maxpairs.lcp,
                                 &TOP_ESA_BOTTOMUP_maxpairs.info,
                                 lastintervallcp,
                                 lastintervalrb - lastintervallb + 1,
                                 NULL,
                                 bustate,
                                 err) != 0)
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
                            &TOP_ESA_BOTTOMUP_maxpairs.info,
#ifdef GT_ESA_BOTTOM_UP_SEQNUM_RELPOS
                            previousseqnum,
                            previousrelpos,
#else
                            previoussuffix,
#endif
                            bustate,
                            err) != 0)
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
    lastsuftabvalue = bucketofsuffixes[numberofsuffixes-1];
    if (processleafedge_maxpairs(false,
                        TOP_ESA_BOTTOMUP_maxpairs.lcp,
                        &TOP_ESA_BOTTOMUP_maxpairs.info,
#ifdef GT_ESA_BOTTOM_UP_SEQNUM_RELPOS
                        gt_seqnumrelpos_decode_seqnum(snrp,lastsuftabvalue),
                        gt_seqnumrelpos_decode_relpos(snrp,lastsuftabvalue),
#else
                        lastsuftabvalue,
#endif
                        bustate,
                        err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_maxpairs.rb = idx;
      if (processlcpinterval_maxpairs(TOP_ESA_BOTTOMUP_maxpairs.lcp,
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
/* vim: ft=c
*/
