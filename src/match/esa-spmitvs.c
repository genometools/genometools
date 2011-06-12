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

#include "core/unused_api.h"
#include "esa-seqread.h"
#include "esa-spmitvs.h"
#include "esa-bottomup.h"

typedef struct  /* global information */
{
  unsigned long unnecessaryleaves, totallength;
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
    return (cc == SEPARATOR) ? true : false;
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
    GtUchar cc = gt_encseq_get_encoded_char(spmitv_state->encseq,
                                            leafnumber + fd,
                                            spmitv_state->readmode);
    if (cc != SEPARATOR &&
        !gt_iswholeleaf(spmitv_state->encseq,spmitv_state->readmode,leafnumber))
    {
      spmitv_state->unnecessaryleaves++;
    }
  }
  return 0;
}

static int processbranchingedge_spmitv(GT_UNUSED bool firstsucc,
                               GT_UNUSED unsigned long fd,
                               GT_UNUSED unsigned long flb,
                               GT_UNUSED GtBUinfo *finfo,
                               GT_UNUSED unsigned long sd,
                               GT_UNUSED unsigned long slb,
                               GT_UNUSED unsigned long srb,
                               GT_UNUSED GtBUinfo *sinfo,
                               GT_UNUSED GtBUstate *bustate,
                               GT_UNUSED GtError *err)
{
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

    state.unnecessaryleaves = 0;
    state.encseq = gt_encseqSequentialsuffixarrayreader(ssar);
    state.readmode = gt_readmodeSequentialsuffixarrayreader(ssar);
    state.totallength = gt_encseq_total_length(state.encseq);
    if (gt_esa_bottomup(ssar, NULL, NULL, processleafedge_spmitv,
                        processbranchingedge_spmitv, (GtBUstate *) &state,
                        err) != 0)
    {
      haserr = true;
    } else
    {
      printf("unnecessaryleaves=%lu (%.2f)\n",
              state.unnecessaryleaves,
              (double) state.unnecessaryleaves/(state.totallength -
                                    gt_encseq_specialcharacters(state.encseq)));
    }
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
