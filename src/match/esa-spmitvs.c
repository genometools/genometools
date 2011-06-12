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
  const GtEncseq *encseq;
} Spmitv_state;

static int processleafedge_spmitv(bool firstsucc,
                          unsigned long fd,
                          unsigned long flb,
                          GT_UNUSED GtBUinfo *info,
                          unsigned long leafnumber,
                          GT_UNUSED GtBUstate *bustate,
                          GT_UNUSED GtError *err)

{
  printf("L %c %lu %lu %lu\n",firstsucc ? '1' : '0',fd,flb,leafnumber);
  return 0;
}

static int processbranchingedge_spmitv(bool firstsucc,
                               unsigned long fd,
                               unsigned long flb,
                               GT_UNUSED GtBUinfo *finfo,
                               unsigned long sd,
                               unsigned long slb,
                               GT_UNUSED unsigned long srb,
                               GT_UNUSED GtBUinfo *sinfo,
                               GT_UNUSED GtBUstate *bustate,
                               GT_UNUSED GtError *err)
{
  printf("B %c %lu %lu %lu %lu\n",firstsucc ? '1' : '0',fd,flb,sd,slb);
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
    if (gt_esa_bottomup(ssar, NULL, NULL, processleafedge_spmitv,
                        processbranchingedge_spmitv, NULL, err) != 0)
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
