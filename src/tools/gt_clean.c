/*
  Copyright (c) 2007-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <string.h>
#include "core/bioseq.h"
#include "core/md5_tab.h"
#include "core/option_api.h"
#include "core/versionfunc.h"
#include "core/xposix.h"
#include "tools/gt_clean.h"

static GtOPrval parse_options(int *parsed_args, int argc, const char **argv,
                              GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("",
                            "Remove all files in the current directory which "
                            "are automatically created by gt.");
  gt_option_parser_set_max_args(op, 0);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
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
  gt_xglob(gt_str_get(path), GLOB_NOCHECK, NULL, &g);

  /* remove found files */
  if (g.gl_pathc) {
    files_to_remove = g.gl_pathv;
    if (strcmp(*files_to_remove, gt_str_get(path))) {
      while (*files_to_remove) {
        gt_xunlink(*files_to_remove);
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
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }
  gt_assert(parsed_args == 1);

  remove_pattern_in_current_dir(GT_ENCSEQFILESUFFIX);
  remove_pattern_in_current_dir(GT_SSPTABFILESUFFIX);
  remove_pattern_in_current_dir(GT_DESTABFILESUFFIX);
  remove_pattern_in_current_dir(GT_SDSTABFILESUFFIX);
  remove_pattern_in_current_dir(GT_OISTABFILESUFFIX);
  remove_pattern_in_current_dir(GT_MD5TABFILESUFFIX);

  return 0;
}
