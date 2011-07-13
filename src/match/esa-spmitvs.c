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
#include "core/mathsupport.h"
#include "esa-seqread.h"
#include "esa-spmitvs.h"
#include "esa-bottomup.h"

typedef struct
{
  unsigned long wholeleaf, wholeleafwidth, nowholeleaf, nowholeleafwidth;
} Lcpintervalcount;

typedef struct  /* global information */
{
  unsigned long unnecessaryleaves,
                totallength,
                currentleafindex,
                lastwholeleaf,
                maxlen;
  unsigned int prefixlength;
  Lcpintervalcount *wholeleafcount;
  const GtEncseq *encseq;
  GtReadmode readmode;
} Spmitv_state;

static bool gt_iswholeleaf(const GtEncseq *encseq,GtReadmode readmode,
                           unsigned long leafnumber)
{
  return (leafnumber > 0)
    ? gt_encseq_position_is_separator(encseq,leafnumber - 1,readmode)
    : true;
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

  if (gt_iswholeleaf(spmitv_state->encseq,spmitv_state->readmode,leafnumber))
  {
    gt_assert(spmitv_state->currentleafindex != spmitv_state->totallength);
    spmitv_state->lastwholeleaf = spmitv_state->currentleafindex;
  } else
  {
    if (leafnumber + fd < spmitv_state->totallength &&
        !gt_encseq_position_is_separator(spmitv_state->encseq,
                                         leafnumber + fd,
                                         spmitv_state->readmode))
    {
      spmitv_state->unnecessaryleaves++;
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
  unsigned long idx;

  for (idx=fd+1; idx<sd; idx++)
  {
    gt_assert(idx <= spmitv_state->maxlen);
    if (spmitv_state->lastwholeleaf != spmitv_state->totallength &&
        spmitv_state->lastwholeleaf >= slb)
    {
      gt_assert(spmitv_state->lastwholeleaf <= srb);
      spmitv_state->wholeleafcount[idx].wholeleaf++;
      spmitv_state->wholeleafcount[idx].wholeleafwidth += (srb - slb + 1);
    } else
    {
      spmitv_state->wholeleafcount[idx].nowholeleaf++;
      spmitv_state->wholeleafcount[idx].nowholeleafwidth += (srb - slb + 1);
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

  if (spmitv_state->lastwholeleaf != spmitv_state->totallength &&
      spmitv_state->lastwholeleaf >= lb)
  {
    gt_assert(lcp <= (unsigned long) spmitv_state->maxlen);
    gt_assert(spmitv_state->lastwholeleaf <= rb);
    spmitv_state->wholeleafcount[lcp].wholeleaf++;
    spmitv_state->wholeleafcount[lcp].wholeleafwidth += (rb - lb + 1);
  } else
  {
    spmitv_state->wholeleafcount[lcp].nowholeleaf++;
    spmitv_state->wholeleafcount[lcp].nowholeleafwidth += (rb - lb + 1);
  }
  return 0;
}

int gt_process_spmitv(const char *inputindex, GtLogger *logger, GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(inputindex,
                                                   SARR_LCPTAB |
                                                   SARR_SUFTAB |
                                                   SARR_ESQTAB,
                                                   SEQ_scan,
                                                   logger,
                                                   err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    Spmitv_state state;
    unsigned int numofchars;
    unsigned long nonspecials;
    GtCodetype numofallcodes;

    state.encseq = gt_encseqSequentialsuffixarrayreader(ssar);
    state.readmode = gt_readmodeSequentialsuffixarrayreader(ssar);
    state.maxlen = gt_encseq_max_seq_length(state.encseq);
    state.unnecessaryleaves = 0;
    state.totallength = gt_encseq_total_length(state.encseq);
    state.currentleafindex = 0;
    state.lastwholeleaf = state.totallength; /* undefined */
    state.prefixlength = gt_Sequentialsuffixarrayreader_prefixlength(ssar);
    state.wholeleafcount = gt_malloc(sizeof (*state.wholeleafcount) *
                                     (state.maxlen+1));
    memset(state.wholeleafcount,0,
           sizeof (*state.wholeleafcount) * (state.maxlen+1));
    numofchars = gt_encseq_alphabetnumofchars(state.encseq);
    numofallcodes = gt_power_for_small_exponents(numofchars,state.prefixlength);
    nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
    if (gt_esa_bottomup(ssar, NULL, NULL, processleafedge_spmitv,
                        processbranchingedge_spmitv,
                        processlcpinterval_spmitv,
                        (GtBUstate *) &state,
                        err) != 0)
    {
      haserr = true;
    } else
    {
      unsigned long idx;

      printf("unnecessaryleaves=%lu (%.2f)\n",
              state.unnecessaryleaves,
              (double) state.unnecessaryleaves/nonspecials);
      for (idx = 0; idx<= state.maxlen; idx++)
      {
        if (state.wholeleafcount[idx].wholeleaf != 0 ||
            state.wholeleafcount[idx].nowholeleaf != 0)
        {
          printf("wholeleaf[%lu]:num=%lu (%.2f), ",idx,
                 state.wholeleafcount[idx].wholeleaf,
                 (double) state.wholeleafcount[idx].wholeleaf/
                          (state.wholeleafcount[idx].wholeleaf+
                           state.wholeleafcount[idx].nowholeleaf));
          printf("width=%lu (%.2f)\n",
                  state.wholeleafcount[idx].wholeleafwidth,
                  (double) state.wholeleafcount[idx].wholeleafwidth/
                           state.totallength);
        }
      }
    }
    gt_free(state.wholeleafcount);
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
