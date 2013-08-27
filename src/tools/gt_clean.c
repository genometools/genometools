/*
  Copyright (c) 2007-2011, 2013 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008       Center for Bioinformatics, University of Hamburg

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
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/xposix.h"
#include "tools/gt_clean.h"

static GtOptionParser* gt_clean_option_parser_new(GT_UNUSED void
                                                  *tool_arguments)
{
  GtOptionParser *op;
  op = gt_option_parser_new("",
                            "Remove all files in the current directory which "
                            "are automatically created by gt.");
  gt_option_parser_set_max_args(op, 0);
  return op;
}

#ifndef _WIN32
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
#endif

int gt_clean_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
                    GT_UNUSED int parsed_args, GT_UNUSED void *tool_arguments,
                    GT_UNUSED GtError *err)
{
  gt_error_check(err);

#ifndef _WIN32
  remove_pattern_in_current_dir(GT_ENCSEQFILESUFFIX);
  remove_pattern_in_current_dir(GT_SSPTABFILESUFFIX);
  remove_pattern_in_current_dir(GT_DESTABFILESUFFIX);
  remove_pattern_in_current_dir(GT_SDSTABFILESUFFIX);
  remove_pattern_in_current_dir(GT_OISTABFILESUFFIX);
  remove_pattern_in_current_dir(GT_MD5TABFILESUFFIX);
#else
  /* XXX */
  gt_error_set(err, "gt_clean_runner() not implemented");
  return -1;
#endif

  return 0;
}

GtTool* gt_clean(void)
{
  return gt_tool_new(NULL,
                     NULL,
                     gt_clean_option_parser_new,
                     NULL,
                     gt_clean_runner);
}
