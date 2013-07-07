/*
  Copyright (c) 2010 Gordon Gremme <gordon@gremme.org>

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
#include "gth/bssm_param_rmsd.h"
#include "gth/gt_gthbssmrmsd.h"

static GtOptionParser* gt_gthbssmrmsd_option_parser_new(GT_UNUSED
                                                        void *tool_arguments)
{
  GtOptionParser *op = gt_option_parser_new("BSSM_file_1 BSSM_file_2", "Show "
                                            "RMSDs between given BSSM files.");
  gt_option_parser_set_mail_address(op, "<gordon@gremme.org>");
  gt_option_parser_set_min_max_args(op, 2, 2);
  return op;
}

static int gt_gthbssmrmsd_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args,
                                 GT_UNUSED void *tool_arguments, GtError *err)
{
  gt_error_check(err);
  return gth_bssm_param_rmsd_show(argv[parsed_args], argv[parsed_args+1], err);
}

GtTool* gt_gthbssmrmsd(void)
{
  return gt_tool_new(NULL,
                     NULL,
                     gt_gthbssmrmsd_option_parser_new,
                     NULL,
                     gt_gthbssmrmsd_runner);
}
