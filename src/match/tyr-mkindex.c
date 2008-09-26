/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/str.h"
#include "core/unused_api.h"
#include "esa-seqread.h"
#include "verbose-def.h"
#include "spacedef.h"
#include "tyr-mkindex.h"

struct Dfsinfo /* information stored for each node of the lcp interval tree */
{
  Seqpos leftmostleaf;
};

struct Dfsstate /* global information */
{
  unsigned long searchlength;
};

#include "esa-dfs.h"

static Dfsinfo *allocateDfsinfo(GT_UNUSED Dfsstate *state)
{
  Dfsinfo *dfsinfo;

  ALLOCASSIGNSPACE(dfsinfo,NULL,Dfsinfo,1);
  return dfsinfo;
}

static void freeDfsinfo(Dfsinfo *dfsinfo, GT_UNUSED Dfsstate *state)
{
  FREESPACE(dfsinfo);
}

static int processleafedge(GT_UNUSED bool firstsucc,
                           GT_UNUSED Seqpos fatherdepth,
                           GT_UNUSED Dfsinfo *father,
                           GT_UNUSED Seqpos leafnumber,
                           GT_UNUSED Dfsstate *state,
                           GT_UNUSED GtError *err)
{
  return 0;
}

static void assignleftmostleaf(Dfsinfo *dfsinfo,Seqpos leftmostleaf,
                               GT_UNUSED Dfsstate *dfsstate)
{
  dfsinfo->leftmostleaf = leftmostleaf;
}

static int enumeratelcpintervals(Sequentialsuffixarrayreader *ssar,
                                 unsigned long searchlength,
                                 Verboseinfo *verboseinfo,
                                 GtError *err)
{
  Dfsstate state;
  bool haserr = false;

  state.searchlength = searchlength;

  if (depthfirstesa(ssar,
                    allocateDfsinfo,
                    freeDfsinfo,
                    processleafedge,
                    NULL,
                    NULL,
                    assignleftmostleaf,
                    NULL,
                    &state,
                    verboseinfo,
                    err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}

int merstatistics(const GtStr *str_inputindex,
                  unsigned long searchlength,
                  GT_UNUSED unsigned long minocc,
                  GT_UNUSED unsigned long maxocc,
                  GT_UNUSED const GtStr *str_storeindex,
                  GT_UNUSED bool storecounts,
                  bool scanfile,
                  bool verbose,
                  GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = newSequentialsuffixarrayreaderfromfile(str_inputindex,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB,
                                                scanfile
                                                  ? SEQ_scan : SEQ_mappedboth,
                                                err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    Verboseinfo *verboseinfo = newverboseinfo(verbose);
    if (enumeratelcpintervals(ssar,
                              searchlength,
                              verboseinfo,
                              err) != 0)
    {
      haserr = true;
    }
    freeverboseinfo(&verboseinfo);
  }
  if (ssar != NULL)
  {
    freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
