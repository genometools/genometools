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
#include "core/logger.h"
#include "core/ma_api.h"
#include "esa-seqread.h"
#include "tyr-occratio.h"

typedef struct /* information stored for each node of the lcp interval tree */
{
  unsigned long leftmostleaf,
                rightmostleaf,
                lcptabrightmostleafplus1;
} OccDfsinfo;

typedef struct /* global information */
{
  const GtEncseq *encseq;
  GtReadmode readmode;
  unsigned long totallength;
  unsigned long minmersize,
                maxmersize;
  GtArrayuint64_t *uniquedistribution,
                  *nonuniquedistribution,
                  *nonuniquemultidistribution;
} OccDfsstate;

#include "esa-dfs.h"

static Dfsinfo* occ_allocateDfsinfo(GT_UNUSED Dfsstate *state)
{
  OccDfsinfo *dfsinfo;

  dfsinfo = gt_malloc(sizeof *dfsinfo);
  return (Dfsinfo *) dfsinfo;
}

static void occ_freeDfsinfo(Dfsinfo *adfsinfo, GT_UNUSED Dfsstate *state)
{
  OccDfsinfo *dfsinfo = (OccDfsinfo*) adfsinfo;
  gt_free(dfsinfo);
}

static void adddistributionuint64_t(GtArrayuint64_t *occdistribution,
                                    unsigned long countocc,
                                    unsigned long value)
{
  if (countocc >= occdistribution->allocateduint64_t)
  {
    const unsigned long addamount = 128UL;
    unsigned long idx;

    occdistribution->spaceuint64_t
      = gt_realloc(occdistribution->spaceuint64_t,
                   sizeof *occdistribution->spaceuint64_t *
                   (countocc+addamount));
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
                                const GtEncseq *encseq,
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
         ISNOTSPECIAL(gt_encseq_get_encoded_char(encseq,pos,readmode));
         pos++, ulen++)
    {
      if (ulen >= (unsigned long) minmersize)
      {
        adddistributionuint64_t(distribution,(unsigned long) ulen,1UL);
      }
    }
  }
}

static int occ_processleafedge(GT_UNUSED bool firstsucc,
                           unsigned long fatherdepth,
                           GT_UNUSED Dfsinfo *father,
                           unsigned long leafnumber,
                           Dfsstate *astate,
                           GT_UNUSED GtError *err)
{
  OccDfsstate *state = (OccDfsstate*) astate;
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

static int occ_processcompletenode(unsigned long nodeptrdepth,
                               Dfsinfo *anodeptr,
                               unsigned long nodeptrminusonedepth,
                               Dfsstate *astate,
                               GT_UNUSED GtError *err)
{
  unsigned long fatherdepth;
  unsigned long startlength, endlength;
  OccDfsinfo *nodeptr = (OccDfsinfo*) anodeptr;
  OccDfsstate *state = (OccDfsstate*) astate;

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

static void occ_assignleftmostleaf(Dfsinfo *adfsinfo,unsigned long leftmostleaf,
                                   GT_UNUSED Dfsstate *dfsstate)
{
  OccDfsinfo *dfsinfo = (OccDfsinfo*) adfsinfo;
  dfsinfo->leftmostleaf = leftmostleaf;
}

static void occ_assignrightmostleaf(Dfsinfo *adfsinfo,
                                    unsigned long currentindex,
                                    GT_UNUSED unsigned long previoussuffix,
                                    unsigned long currentlcp,
                                    GT_UNUSED Dfsstate *dfsstate)
{
  OccDfsinfo *dfsinfo = (OccDfsinfo*) adfsinfo;
  dfsinfo->rightmostleaf = currentindex;
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
  OccDfsstate *state;
  bool haserr = false;

  gt_error_check(err);
  state = gt_malloc(sizeof (*state));
  state->encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  state->readmode = gt_readmodeSequentialsuffixarrayreader(ssar);
  state->totallength = gt_encseq_total_length(state->encseq);
  state->minmersize = minmersize;
  state->maxmersize = maxmersize;
  state->uniquedistribution = uniquedistribution;
  state->nonuniquedistribution = nonuniquedistribution;
  state->nonuniquemultidistribution = nonuniquemultidistribution;
  if (gt_depthfirstesa(ssar,
                    occ_allocateDfsinfo,
                    occ_freeDfsinfo,
                    occ_processleafedge,
                    NULL,
                    occ_processcompletenode,
                    occ_assignleftmostleaf,
                    occ_assignrightmostleaf,
                    (Dfsstate*) state,
                    logger,
                    err) != 0)
  {
    haserr = true;
  }
  gt_free(state);
  return haserr ? -1 : 0;
}

int gt_tyr_occratio_func(const char *inputindex,
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
  ssar = gt_newSequentialsuffixarrayreaderfromfile(inputindex,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB,
                                                scanfile ? SEQ_scan
                                                         : SEQ_mappedboth,
                                                logger,
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
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
