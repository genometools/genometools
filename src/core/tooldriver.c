/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include <stdlib.h>
#include "core/init.h"
#include "core/error.h"
#include "core/tooldriver.h"

int gt_tooldriver(GtToolFunc tool_func, int argc, char *argv[])
{
  GtError *err;
  int had_err;
  gt_lib_init();
  err = gt_error_new();
  gt_error_set_progname(err, argv[0]);
  had_err = tool_func(argc, (const char**) argv, err);
  if (gt_error_is_set(err)) {
    fprintf(stderr, "%s: error: %s\n", gt_error_get_progname(err),
            gt_error_get(err));
    gt_assert(had_err);
  }
  gt_error_delete(err);
  if (gt_lib_clean())
    return 2; /* programmer error */
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

int gt_toolobjdriver(GtToolConstructor tool_constructor, int argc, char *argv[])
{
  GtTool *tool;
  GtError *err;
  int had_err;
  gt_lib_init();
  err = gt_error_new();
  gt_error_set_progname(err, argv[0]);
  tool = tool_constructor();
  had_err = gt_tool_run(tool, argc, (const char**) argv, err);
  gt_tool_delete(tool);
  if (gt_error_is_set(err)) {
    fprintf(stderr, "%s: error: %s\n", gt_error_get_progname(err),
            gt_error_get(err));
    gt_assert(had_err);
  }
  gt_error_delete(err);
  if (gt_lib_clean())
    return 2; /* programmer error */
  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
