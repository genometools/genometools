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

#include "core/option_api.h"
#include "core/versionfunc.h"
#include "gth/bssm_param.h"
#include "gth/gthspeciestab.h" /* XXX */
#include "gth/gt_gthmkbssmfiles.h"

/* The number of the last species defined in header file */
#define LASTSPECIESNUM 9

static GtOPrval gthmkbssmfiles_parse_options(int *parsed_args, int argc,
                                             const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("output_path", "Write hard coded BSSM files to "
                         "output_path.");
  gt_option_parser_set_min_max_args(op, 1, 1);
  gt_option_parser_set_mail_address(op, "<gremme@zbh.uni-hamburg.de>");
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_gthmkbssmfiles(int argc, const char **argv, GtError *err)
{
  unsigned long i;
  GtStr *filename;
  int parsed_args, had_err = 0;

  /* option parsing */
  switch (gthmkbssmfiles_parse_options(&parsed_args, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }

  gt_assert(parsed_args + 1 == argc);
  filename = gt_str_new();

  for (i = 0; !had_err && i <= LASTSPECIESNUM; i++) {
    GthBSSMParam *bssm_param;
    gt_str_append_cstr(filename, argv[parsed_args]);
    gt_str_append_char(filename, '/');
    gt_str_append_cstr(filename, speciestab[i]);

    /* for files which are obsolete due to new model files produced by
       gthbssmbuild add an .old after the species name */
    if (i >= 8)
      gt_str_append_cstr(filename, ".old");

    gt_str_append_char(filename, '.');
    gt_str_append_cstr(filename, BSSMFILEENDING);

    if (!(bssm_param = gth_bssm_param_extract(i, err)))
      had_err = -1;
    if (!had_err)
      had_err = gth_bssm_param_save(bssm_param, gt_str_get(filename), err);
    gth_bssm_param_delete(bssm_param);

    /* resetting filename */
    gt_str_reset(filename);
  }

  gt_str_delete(filename);

  return had_err;
}
