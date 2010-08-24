/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "core/unused_api.h"
#include "core/logger.h"
#include "esa-seqread.h"

typedef struct /* information stored for each node of the lcp interval tree */
{
  unsigned long *filenumdist;
} ShulengthdistDfsinfo;

typedef struct  /* global information */
{
  unsigned long numofdbfiles;
  const GtEncseq *encseq;
} Shulengthdiststate;

#include "esa-dfs.h"

static Dfsinfo *allocateDfsinfo(Dfsstate *astate)
{
  ShulengthdistDfsinfo *dfsinfo;
  Shulengthdiststate *state = (Shulengthdiststate*) astate;

  dfsinfo = gt_malloc(sizeof(*dfsinfo));
  dfsinfo->filenumdist 
    = gt_malloc(sizeof(*dfsinfo->filenumdist) * state->numofdbfiles);
  return (Dfsinfo*) dfsinfo;
}

static void freeDfsinfo(Dfsinfo *adfsinfo, GT_UNUSED Dfsstate *state)
{
  ShulengthdistDfsinfo *dfsinfo = (ShulengthdistDfsinfo*) adfsinfo;;
  gt_free(dfsinfo->filenumdist);
  gt_free(dfsinfo);
}

static int processleafedge(GT_UNUSED bool firstsucc,
                           unsigned long fatherdepth,
                           GT_UNUSED Dfsinfo *afather,
                           GT_UNUSED unsigned long leafnumber,
                           GT_UNUSED Dfsstate *astate,
                           GT_UNUSED GtError *err)
{
  /*
  Shulengthdiststate *state = (Shulengthdiststate*) astate;
  ShulengthdistDfsinfo *father = (ShulengthdistDfsinfo*) afather;
  */

#ifdef SKDEBUG
  printf("processleafedge %lu firstsucc=%s, "
         " depth(father)= %lu\n",
         leafnumber,
         firstsucc ? "true" : "false",
         fatherdepth);
#endif
  if (fatherdepth == 0)
  {
    return 0;
  }
  return 0;
}

static int processbranchedge(GT_UNUSED bool firstsucc,
                             unsigned long fatherdepth,
                             GT_UNUSED Dfsinfo *afather,
                             GT_UNUSED Dfsinfo *ason,
                             GT_UNUSED Dfsstate *astate,
                             GT_UNUSED GtError *err)
{
  /*
  Shulengthdiststate *state = (Shulengthdiststate*) astate;
  ShulengthdistDfsinfo *son = (ShulengthdistDfsinfo*) ason;
  ShulengthdistDfsinfo *father = (ShulengthdistDfsinfo*) afather;
  */

#ifdef SKDEBUG
  printf("processbranchedge firstsucc=%s, depth(father)=%lu\n",
         firstsucc ? "true" : "false",fatherdepth);
#endif
  if (fatherdepth == 0)
  {
    return 0;
  }
  return 0;
}

int gt_sastream2shulengthdist(Sequentialsuffixarrayreader *ssar,
                              const GtEncseq *encseq,
                              GtLogger *logger,
                              GtError *err)
{
  Shulengthdiststate *state;
  bool haserr = false;

  state = gt_malloc(sizeof(*state));
  state->numofdbfiles = gt_encseq_num_of_files(encseq);
  state->encseq = encseq;

  if (gt_depthfirstesa(ssar,
                       allocateDfsinfo,
                       freeDfsinfo,
                       processleafedge,
                       processbranchedge,
                       NULL,
                       NULL,
                       NULL,
                       (Dfsstate*) state,
                       logger,
                       err) != 0)
  {
    haserr = true;
  }
  gt_free(state);
  return haserr ? -1 : 0;
}
