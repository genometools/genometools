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

#include <math.h>
#include "core/unused_api.h"
#include "esa-seqread.h"
#include "esa-spmitvs.h"
#include "bcktab.h"
#include "esa-bottomup.h"

typedef struct  /* global information */
{
  unsigned long unnecessaryleaves,
                totallength,
                currentleafindex,
                lastwholeleaf,
                bucketwithwholeleaf,
                bucketwithwholeleafwidth,
                allbuckets;
  unsigned int prefixlength;
  const GtEncseq *encseq;
  GtReadmode readmode;
} Spmitv_state;

static bool gt_iswholeleaf(const GtEncseq *encseq,GtReadmode readmode,
                           unsigned long leafnumber)
{
  if (leafnumber > 0)
  {
    GtUchar cc = gt_encseq_get_encoded_char(encseq,
                                            leafnumber - 1,
                                            readmode);
    return (cc == (GtUchar) SEPARATOR) ? true : false;
  }
  return true;
}

static int processleafedge_spmitv(GT_UNUSED bool firstsucc,
                                  unsigned long fd,
                                  GT_UNUSED unsigned long flb,
                                  GT_UNUSED GtBUinfo *info,
                                  unsigned long leafnumber,
                                  GtBUstate *bustate,
                                  GT_UNUSED GtError *err)

{
  Spmitv_state *spmitv_state = (Spmitv_state *) bustate;

  if (leafnumber + fd < spmitv_state->totallength)
  {
    if (gt_iswholeleaf(spmitv_state->encseq,
                       spmitv_state->readmode,leafnumber))
    {
      gt_assert(spmitv_state->currentleafindex != spmitv_state->totallength);
      spmitv_state->lastwholeleaf = spmitv_state->currentleafindex;
    } else
    {
      GtUchar cc = gt_encseq_get_encoded_char(spmitv_state->encseq,
                                              leafnumber + fd,
                                              spmitv_state->readmode);
      if (cc != (GtUchar) SEPARATOR)
      {
        spmitv_state->unnecessaryleaves++;
      }
    }
  }
  spmitv_state->currentleafindex++;
  return 0;
}

static int processbranchingedge_spmitv(GT_UNUSED bool firstsucc,
                                       unsigned long fd,
                                       GT_UNUSED unsigned long flb,
                                       GT_UNUSED GtBUinfo *finfo,
                                       unsigned long sd,
                                       unsigned long slb,
                                       unsigned long srb,
                                       GT_UNUSED GtBUinfo *sinfo,
                                       GtBUstate *bustate,
                                       GT_UNUSED GtError *err)
{
  Spmitv_state *spmitv_state = (Spmitv_state *) bustate;

  if (fd < (unsigned long) spmitv_state->prefixlength &&
      sd > (unsigned long) spmitv_state->prefixlength)
  {
    spmitv_state->allbuckets++;
    if (spmitv_state->lastwholeleaf != spmitv_state->totallength &&
        spmitv_state->lastwholeleaf >= slb)
    {
      gt_assert(spmitv_state->lastwholeleaf <= srb);
      spmitv_state->bucketwithwholeleaf++;
      spmitv_state->bucketwithwholeleafwidth += (srb - slb + 1);
    }
  }
  return 0;
}

static int processlcpinterval_spmitv(unsigned long lcp,
                                     unsigned long lb,
                                     unsigned long rb,
                                     GT_UNUSED GtBUinfo *info,
                                     GtBUstate *bustate,
                                     GT_UNUSED GtError *err)
{
  Spmitv_state *spmitv_state = (Spmitv_state *) bustate;

  if (lcp == (unsigned long) spmitv_state->prefixlength)
  {
    spmitv_state->allbuckets++;
    if (spmitv_state->lastwholeleaf != spmitv_state->totallength &&
        spmitv_state->lastwholeleaf >= lb)
    {
      gt_assert(spmitv_state->lastwholeleaf <= rb);
      spmitv_state->bucketwithwholeleaf++;
      spmitv_state->bucketwithwholeleafwidth += (rb - lb + 1);
    }
  }
  return 0;
}

int gt_process_spmitv(const char *inputindex, GT_UNUSED unsigned int minlen,
                      GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(inputindex,
                                                   SARR_LCPTAB |
                                                   SARR_SUFTAB |
                                                   SARR_BCKTAB |
                                                   SARR_ESQTAB,
                                                   SEQ_scan,
                                                   err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    Spmitv_state state;
    unsigned int numofchars;
    unsigned long emptybuckets, nonspecials;
    Bcktab *bcktab;
    GtCodetype numofallcodes;

    state.unnecessaryleaves = 0;
    state.encseq = gt_encseqSequentialsuffixarrayreader(ssar);
    state.readmode = gt_readmodeSequentialsuffixarrayreader(ssar);
    state.totallength = gt_encseq_total_length(state.encseq);
    state.currentleafindex = 0;
    state.lastwholeleaf = state.totallength; /* undefined */
    state.allbuckets = 0;
    state.prefixlength = gt_Sequentialsuffixarrayreader_prefixlength(ssar);
    state.bucketwithwholeleaf = 0;
    state.bucketwithwholeleafwidth = 0;
    numofchars = gt_encseq_alphabetnumofchars(state.encseq);
    bcktab = gt_Sequentialsuffixarrayreader_bcktab(ssar);
    numofallcodes = gt_bcktab_numofallcodes(bcktab);
    nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);

    gt_determinemaxbucketsize(bcktab,
                              0,
                              numofallcodes-1,
                              nonspecials,
                              numofchars);
    emptybuckets = gt_bcktab_emptybuckets(bcktab);
    printf("emptybuckets=%lu (%.2f)\n",emptybuckets,
                                       (double) emptybuckets/numofallcodes);
    if (gt_esa_bottomup(ssar, NULL, NULL, processleafedge_spmitv,
                        processbranchingedge_spmitv,
                        processlcpinterval_spmitv,
                        (GtBUstate *) &state,
                        err) != 0)
    {
      haserr = true;
    } else
    {
      printf("unnecessaryleaves=%lu (%.2f)\n",
              state.unnecessaryleaves,
              (double) state.unnecessaryleaves/(state.totallength -
                                    gt_encseq_specialcharacters(state.encseq)));
      printf("allbuckets=%lu (%.2f)\n",state.allbuckets,
                                       (double) state.allbuckets/numofallcodes);
      printf("bucketwithwholeleaf=%lu (%.2f)\n",state.bucketwithwholeleaf,
              (double) state.bucketwithwholeleaf/(numofallcodes-emptybuckets));
      printf("bucketwithwholeleafwidth=%lu (%.2f)\n",
              state.bucketwithwholeleafwidth,
              (double) state.bucketwithwholeleafwidth/nonspecials);
    }
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
