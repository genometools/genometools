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
#include "core/logger.h"
#include "spacedef.h"
#include "tyr-occratio.h"

struct Dfsinfo /* information stored for each node of the lcp interval tree */
{
  unsigned long leftmostleaf,
         rightmostleaf,
         suftabrightmostleaf,
         lcptabrightmostleafplus1;
};

struct Dfsstate /* global information */
{
  const GtEncodedsequence *encseq;
  GtReadmode readmode;
  unsigned long totallength;
  unsigned long minmersize,
                maxmersize;
  GtArrayuint64_t *uniquedistribution,
                  *nonuniquedistribution,
                  *nonuniquemultidistribution;
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

static void adddistributionuint64_t(GtArrayuint64_t *occdistribution,
                                    unsigned long countocc,
                                    unsigned long value)
{
  if (countocc >= occdistribution->allocateduint64_t)
  {
    const unsigned long addamount = 128UL;
    unsigned long idx;

    ALLOCASSIGNSPACE(occdistribution->spaceuint64_t,
                     occdistribution->spaceuint64_t,
                     uint64_t,countocc+addamount);
    for (idx=occdistribution->allocateduint64_t;
         idx<countocc+addamount; idx++)
    {
      occdistribution->spaceuint64_t[idx] = 0;
    }
    occdistribution->allocateduint64_t = countocc+addamount;
  }
  if (countocc + 1 > occdistribution->nextfreeuint64_t)
  {
    occdistribution->nextfreeuint64_t = countocc+1;
  }
  occdistribution->spaceuint64_t[countocc] += value;
}

static void iteritvdistribution(GtArrayuint64_t *distribution,
                                const GtEncodedsequence *encseq,
                                GtReadmode readmode,
                                unsigned long totallength,
                                unsigned long minmersize,
                                unsigned long maxmersize,
                                unsigned long length,
                                unsigned long startpos)
{

  if (length <= (unsigned long) maxmersize)
  {
    unsigned long ulen, pos;

    for (ulen = length,
         pos = startpos + length - 1;
         ulen <= (unsigned long) maxmersize &&
         pos < totallength &&
         ISNOTSPECIAL(gt_encodedsequence_getencodedchar(encseq,pos,readmode));
         pos++, ulen++)
    {
      if (ulen >= (unsigned long) minmersize)
      {
        adddistributionuint64_t(distribution,(unsigned long) ulen,1UL);
      }
    }
  }
}

static int processleafedge(GT_UNUSED bool firstsucc,
                           unsigned long fatherdepth,
                           GT_UNUSED Dfsinfo *father,
                           unsigned long leafnumber,
                           Dfsstate *state,
                           GT_UNUSED GtError *err)
{
  iteritvdistribution(state->uniquedistribution,
                      state->encseq,
                      state->readmode,
                      state->totallength,
                      state->minmersize,
                      state->maxmersize,
                      fatherdepth+1,
                      leafnumber);
  return 0;
}

static int processcompletenode(unsigned long nodeptrdepth,
                               Dfsinfo *nodeptr,
                               unsigned long nodeptrminusonedepth,
                               Dfsstate *state,
                               GT_UNUSED GtError *err)
{
  unsigned long fatherdepth;
  unsigned long startlength, endlength;

  fatherdepth = nodeptr->lcptabrightmostleafplus1;
  if (fatherdepth < nodeptrminusonedepth)
  {
    fatherdepth = nodeptrminusonedepth;
  }
  startlength = (unsigned long) (fatherdepth + 1);
  if (startlength < state->minmersize)
  {
    startlength = state->minmersize;
  }
  endlength = (unsigned long) nodeptrdepth;
  if (endlength > state->maxmersize)
  {
    endlength = state->maxmersize;
  }
  if (startlength <= endlength)
  {
    unsigned long lenval;
    unsigned long occcount = nodeptr->rightmostleaf - nodeptr->leftmostleaf + 1;

    for (lenval = startlength; lenval <= endlength; lenval++)
    {
      adddistributionuint64_t(state->nonuniquemultidistribution,
                              lenval,
                              (unsigned long) occcount);
      adddistributionuint64_t(state->nonuniquedistribution,
                              lenval,
                              1UL);
    }
  }
  return 0;
}

static void assignleftmostleaf(Dfsinfo *dfsinfo,unsigned long leftmostleaf,
                               GT_UNUSED Dfsstate *dfsstate)
{
  dfsinfo->leftmostleaf = leftmostleaf;
}

static void assignrightmostleaf(Dfsinfo *dfsinfo,unsigned long currentindex,
                                unsigned long previoussuffix,
                                unsigned long currentlcp,
                                GT_UNUSED Dfsstate *dfsstate)
{
  dfsinfo->rightmostleaf = currentindex;
  dfsinfo->suftabrightmostleaf = previoussuffix;
  dfsinfo->lcptabrightmostleafplus1 = currentlcp;
}

static int computeoccurrenceratio(Sequentialsuffixarrayreader *ssar,
                                  unsigned long minmersize,
                                  unsigned long maxmersize,
                                  GtArrayuint64_t *uniquedistribution,
                                  GtArrayuint64_t *nonuniquedistribution,
                                  GtArrayuint64_t *nonuniquemultidistribution,
                                  GtLogger *logger,
                                  GtError *err)
{
  Dfsstate state;
  bool haserr = false;

  gt_error_check(err);
  state.encseq = encseqSequentialsuffixarrayreader(ssar);
  state.readmode = readmodeSequentialsuffixarrayreader(ssar);
  state.totallength = gt_encodedsequence_total_length(state.encseq);
  state.minmersize = minmersize;
  state.maxmersize = maxmersize;
  state.uniquedistribution = uniquedistribution;
  state.nonuniquedistribution = nonuniquedistribution;
  state.nonuniquemultidistribution = nonuniquemultidistribution;
  if (depthfirstesa(ssar,
                    allocateDfsinfo,
                    freeDfsinfo,
                    processleafedge,
                    NULL,
                    processcompletenode,
                    assignleftmostleaf,
                    assignrightmostleaf,
                    &state,
                    logger,
                    err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}

int tyr_occratio(const GtStr *str_inputindex,
                 bool scanfile,
                 unsigned long minmersize,
                 unsigned long maxmersize,
                 GtArrayuint64_t *uniquedistribution,
                 GtArrayuint64_t *nonuniquedistribution,
                 GtArrayuint64_t *nonuniquemultidistribution,
                 GtLogger *logger,
                 GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = newSequentialsuffixarrayreaderfromfile(str_inputindex,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB,
                                                scanfile ? SEQ_scan
                                                         : SEQ_mappedboth,
                                                err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (computeoccurrenceratio(ssar,
                               minmersize,
                               maxmersize,
                               uniquedistribution,
                               nonuniquedistribution,
                               nonuniquemultidistribution,
                               logger,
                               err) != 0)
    {
      haserr = true;
    }
  }
  if (ssar != NULL)
  {
    freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
