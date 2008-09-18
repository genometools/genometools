/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/bioseq.h"
#include "core/option.h"
#include "core/versionfunc.h"
#include "core/xposix.h"
#include "tools/gt_clean.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            GtError *err)
{
  GtOptionParser *op;
  OPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("",
                            "Remove all files in the current directory which "
                            "are automatically created by gt.");
  gt_option_parser_set_max_args(op, 0);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

static void remove_pattern_in_current_dir(const char *pattern)
{
  char **files_to_remove;
  GtStr *path;
  glob_t g;

  path = gt_str_new_cstr("./*");
  gt_str_append_cstr(path, pattern);
  xglob(gt_str_get(path), GLOB_NOCHECK, NULL, &g);

  /* remove found files */
  if (g.gl_pathc) {
    files_to_remove = g.gl_pathv;
    if (strcmp(*files_to_remove, gt_str_get(path))) {
      while (*files_to_remove) {
        xunlink(*files_to_remove);
        files_to_remove++;
      }
    }
  }

  /* free */
  globfree(&g);
  gt_str_delete(path);
}

int gt_clean(int argc, const char **argv, GtError *err)
{
  int parsed_args;
  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  /* remove GT_BIOSEQ_INDEX files */
  remove_pattern_in_current_dir(GT_BIOSEQ_INDEX);

  /* remove GT_BIOSEQ_RAW files */
  remove_pattern_in_current_dir(GT_BIOSEQ_RAW);

  /* remove GT_BIOSEQ_FINGERPRINTS files */
  remove_pattern_in_current_dir(GT_BIOSEQ_FINGERPRINTS);

  return 0;
}
