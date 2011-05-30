/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/option_api.h"
#include "gth/bssm_param.h"
#include "gth/gt_gthbssmfileinfo.h"

static GtOPrval gthbssmfileinfo_parse_options(int *parsed_args, int argc,
                                              const char **argv,
                                              GtShowVersionFunc
                                              gth_version_func, GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("bssm_file", "Show information about the specified "
                         "BSSM file.");
  gt_option_parser_set_min_args(op, 1);
  gt_option_parser_set_mail_address(op, "<gremme@zbh.uni-hamburg.de>");
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gth_version_func,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_gthbssmfileinfo(int argc, const char **argv,
                       GtShowVersionFunc gth_version_func, GtError *err)
{
  GthBSSMParam *bssm_param;
  GtStr *bssm_param_filename;
  int parsed_args, had_err = 0;

  gt_error_check(err);

  /* option parsing */
  switch (gthbssmfileinfo_parse_options(&parsed_args, argc, argv,
                                        gth_version_func, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }

  /* build bssm file name */
  bssm_param_filename = gt_str_new_cstr(argv[parsed_args]);
  gt_str_append_char(bssm_param_filename, '.');
  gt_str_append_cstr(bssm_param_filename, BSSMFILEENDING);

  /* load bssm file */
  if (!(bssm_param = gth_bssm_param_load(gt_str_get(bssm_param_filename), err)))
  had_err = -1;

  /* output bssm file info on stdout */
  if (!had_err)
    gth_bssm_param_show_info(bssm_param, NULL);

  gt_str_delete(bssm_param_filename);
  gth_bssm_param_delete(bssm_param);

  return had_err;
}
