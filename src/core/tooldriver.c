/*
  Copyright (c) 2007-2011 Gordon Gremme <gordon@gremme.org>
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

#include <stdio.h>
#include "core/init_api.h"
#include "core/tool.h"
#include "core/tooldriver_api.h"

int gt_tooldriver(GtToolFunc tool, int argc, char *argv[])
{
  gt_assert(tool && argv);
  GtError *err;
  int had_err;
  gt_lib_init();
  gt_assert(tool && argv);
  err = gt_error_new();
  gt_error_set_progname(err, argv[0]);
  had_err = tool(argc, (const char**) argv, err);
  if (gt_error_is_set(err)) {
    fprintf(stderr, "%s: error: %s\n", gt_error_get_progname(err),
            gt_error_get(err));
    gt_assert(had_err);
  }
  gt_error_delete(err);
  if (gt_lib_clean())
    return GT_EXIT_PROGRAMMING_ERROR; /* programmer error */
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

int gt_toolobjdriver(GtToolConstructor tool_constructor,
                     GtShowVersionFunc version_func, int argc, char *argv[])
{
  gt_assert(tool_constructor && argv);
  GtTool *tool;
  GtError *err;
  int had_err;
  gt_lib_init();
  gt_assert(tool_constructor && argv);
  err = gt_error_new();
  gt_error_set_progname(err, argv[0]);
  tool = tool_constructor();
  if (version_func) {
    GtOptionParser *op = gt_tool_get_option_parser(tool);
    gt_option_parser_set_version_func(op, version_func);
  }
  had_err = gt_tool_run(tool, argc, (const char**) argv, err);
  gt_tool_delete(tool);
  if (gt_error_is_set(err)) {
    fprintf(stderr, "%s: error: %s\n", gt_error_get_progname(err),
            gt_error_get(err));
    gt_assert(had_err);
  }
  gt_error_delete(err);
  if (gt_lib_clean())
    return GT_EXIT_PROGRAMMING_ERROR; /* programmer error */
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
