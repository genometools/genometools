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

#include "core/logger.h"
#include "core/unused_api.h"
#include "esa-seqread.h"
#include "esa-lcpintervals.h"

typedef struct /* information stored for each node of the lcp interval tree */
{
  unsigned long leftmostleaf,
                rightmostleaf;
} Elcpinfo;

typedef struct
{
  unsigned long lcpvalue, lb, rb;
} Lcpinterval;

typedef struct  /* global information */
{
  int (*processlcpinterval)(void *,Lcpinterval *);
  void *processinfo;
} Elcpstate;

#include "esa-dfs.h"

static Dfsinfo *elcp_allocateDfsinfo(GT_UNUSED Dfsstate *astate)
{
  return (Dfsinfo *) gt_malloc(sizeof (Elcpinfo));
}

static void elcp_freeDfsinfo(Dfsinfo *adfsinfo,GT_UNUSED Dfsstate *state)
{
  gt_free((Elcpinfo *) adfsinfo);
}

static int elcp_processcompletenode(
                          unsigned long nodeptrdepth,
                          Dfsinfo *anodeptr,
                          GT_UNUSED unsigned long nodeptrminusonedepth,
                          Dfsstate *astate,
                          GT_UNUSED GtError *err)
{
  Elcpinfo *nodeptr = (Elcpinfo *) anodeptr;
  Elcpstate *state = (Elcpstate *) astate;
  Lcpinterval lcpinterval;

  gt_assert(state != NULL);
  gt_assert(nodeptr != NULL);
  lcpinterval.lcpvalue = nodeptrdepth;
  lcpinterval.lb = nodeptr->leftmostleaf;
  lcpinterval.rb = nodeptr->rightmostleaf;
  if (state->processlcpinterval != NULL)
  {
    if (state->processlcpinterval(state->processinfo,&lcpinterval) != 0)
    {
      return -1;
    }
  }
  return 0;
}

static void elcp_assignleftmostleaf(Dfsinfo *adfsinfo,
                                    unsigned long leftmostleaf,
                                    GT_UNUSED Dfsstate *dfsstate)
{
  ((Elcpinfo *) adfsinfo)->leftmostleaf = leftmostleaf;
}

static void elcp_assignrightmostleaf(Dfsinfo *adfsinfo,
                                     unsigned long currentindex,
                                     GT_UNUSED unsigned long previoussuffix,
                                     GT_UNUSED unsigned long currentlcp,
                                     GT_UNUSED Dfsstate *dfsstate)
{
  ((Elcpinfo *) adfsinfo)->rightmostleaf = currentindex;
}

static int gt_enumlcpvalues(Sequentialsuffixarrayreader *ssar,
                            int (*processlcpinterval)(void *,Lcpinterval *),
                            void *processinfo,
                            GtLogger *logger,
                            GtError *err)
{
  Elcpstate *state;
  bool haserr = false;

  state = gt_malloc(sizeof (*state));
  state->processlcpinterval = processlcpinterval;
  state->processinfo = processinfo;
  if (gt_depthfirstesa(ssar,
                       elcp_allocateDfsinfo,
                       elcp_freeDfsinfo,
                       NULL,
                       NULL,
                       elcp_processcompletenode,
                       elcp_assignleftmostleaf,
                       elcp_assignrightmostleaf,
                       (Dfsstate *) state,
                       logger,
                       err) != 0)
  {
    haserr = true;
  }
  gt_free(state);
  return haserr ? -1 : 0;
}

static int showlcpinterval(GT_UNUSED void *data,Lcpinterval *lcpinterval)
{
  printf("%lu %lu %lu\n",lcpinterval->lcpvalue,lcpinterval->lb,lcpinterval->rb);
  return 0;
}

int gt_runenumlcpvalues(const char *inputindex,
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
  if (!haserr)
  {
    if (gt_enumlcpvalues(ssar, showlcpinterval, NULL, logger, err) != 0)
    {
      haserr = true;
    }
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
