/*
  Copyright (c) 2004-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004      Michael E Sparks <mespar1@iastate.edu>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

/*
  This is simply a driver program for using the echo_bssm function to help
  potentially debug our parameterization routines. Namely, it allows us to
  examine the contents of the parameterization without modifying them.

  I think that this is a handy tool to have around, and potentially something
  useful to share with the end user! (For printing ASCII parameterizations that
  would be used in GeneSeqer, for example.)
*/

#include "core/option_api.h"
#include "gth/bssm_param.h"
#include "gth/gt_gthbssmprint.h"

static GtOPrval gthbssmprint_parse_options(int *parsed_args, int argc,
                                           const char **argv,
                                           GtShowVersionFunc gth_version_func,
                                           GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("bssm_file",
                            "Print BSSM file bssm_file to stdout.");
  gt_option_parser_set_min_max_args(op, 1, 1);
  gt_option_parser_set_mail_address(op, "<gremme@zbh.uni-hamburg.de>");
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gth_version_func,
                               err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_gthbssmprint(int argc, const char **argv,
                    GtShowVersionFunc gth_version_func, GtError *err)
{
  GthBSSMParam *bssm_param; /* stores model parameterization */
  int parsed_args, had_err = 0;

  /* verify command line specification of a parameter file to inspect */
  switch (gthbssmprint_parse_options(&parsed_args, argc, argv, gth_version_func,
                                     err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }
  gt_assert(parsed_args + 1 == argc);

  if (!(bssm_param = gth_bssm_param_load(argv[parsed_args], err)))
    had_err = -1;

  /* print debugging information to stdout */
  if (!had_err)
    gth_bssm_param_echo(bssm_param, stdout);

  gth_bssm_param_delete(bssm_param);

  return had_err;
}
