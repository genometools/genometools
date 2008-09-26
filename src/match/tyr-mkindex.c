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
  Seqpos depth,
         leftmostleaf,
         rightmostleaf,
         suftabrightmostleaf,
         lcptabrightmostleafplus1;
};

typedef int (*Processoccurrencecount)(Seqpos,Seqpos,void *);

struct Dfsstate /* global information */
{
  Seqpos searchlength,
         totallength;
  const Encodedsequence *encseq;
  Readmode readmode;
  Processoccurrencecount processoccurrencecount;
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

static bool containsspecial(GT_UNUSED const Encodedsequence *encseq,
                            GT_UNUSED Seqpos startindex,
                            GT_UNUSED Seqpos len)
{
  return false;
}

static int processleafedge(GT_UNUSED bool firstsucc,
                           GT_UNUSED Seqpos fatherdepth,
                           Dfsinfo *father,
                           Seqpos leafnumber,
                           GT_UNUSED Dfsstate *state,
                           GT_UNUSED GtError *err)
{
  if (father->depth < state->searchlength &&
      leafnumber + state->searchlength <=
      state->totallength &&
      !containsspecial(state->encseq,
                      leafnumber + father->depth,
                      state->searchlength - father->depth))
  {
    if (state->processoccurrencecount((Seqpos) 1,leafnumber,state) != 0)
    {
      return -1;
    }
  }
  return 0;
}

static void assignleftmostleaf(Dfsinfo *dfsinfo,Seqpos leftmostleaf,
                               GT_UNUSED Dfsstate *dfsstate)
{
  dfsinfo->leftmostleaf = leftmostleaf;
}

static void assignrightmostleaf(Dfsinfo *dfsinfo,Seqpos currentindex,
                                Seqpos previoussuffix,Seqpos currentlcp,
                                GT_UNUSED Dfsstate *dfsstate)
{
  dfsinfo->rightmostleaf = currentindex;
  dfsinfo->suftabrightmostleaf = previoussuffix;
  dfsinfo->lcptabrightmostleafplus1 = currentlcp;
}

#define ASSIGNRIGHTMOSTLEAF(STACKELEM,VALUE,SUFVALUE,LCPVALUE)\

static int enumeratelcpintervals(Sequentialsuffixarrayreader *ssar,
                                 unsigned long searchlength,
                                 Verboseinfo *verboseinfo,
                                 GtError *err)
{
  Dfsstate state;
  bool haserr = false;

  state.searchlength = (Seqpos) searchlength;
  state.encseq = encseqSequentialsuffixarrayreader(ssar);
  state.totallength = getencseqtotallength(state.encseq);
  state.readmode = readmodeSequentialsuffixarrayreader(ssar);
  if (depthfirstesa(ssar,
                    allocateDfsinfo,
                    freeDfsinfo,
                    processleafedge,
                    NULL,
                    NULL,
                    assignleftmostleaf,
                    assignrightmostleaf,
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
