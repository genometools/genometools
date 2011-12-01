/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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
#include "core/mathsupport.h"
#include "match/esa_spmitvs_visitor.h"
#include "esa-spmitvs.h"
#include "esa-bottomup.h"

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
    GtESAVisitor *ev;
    unsigned long nonspecials;

    nonspecials = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
    ev = gt_esa_spmitvs_visitor_new(
                              gt_encseqSequentialsuffixarrayreader(ssar),
                              gt_readmodeSequentialsuffixarrayreader(ssar),
                              gt_Sequentialsuffixarrayreader_prefixlength(ssar),
                              err);

    if (gt_esa_bottomup(ssar, ev, err) != 0)
    {
      haserr = true;
    } else
    {
      gt_esa_spmitvs_visitor_print_results((GtESASpmitvsVisitor*) ev,
                                           nonspecials);
    }
    gt_esa_visitor_delete(ev);
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
