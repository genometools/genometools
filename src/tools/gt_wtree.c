/*
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include "core/cstr_array.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/toolbox_api.h"
#include "core/unused_api.h"
#include "tools/gt_wtree.h"
#include "tools/gt_wtree_bench.h"

static void* gt_wtree_arguments_new(void)
{
  GtToolbox *wtree_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(wtree_toolbox, "benchmark", gt_wtree_bench());
  return wtree_toolbox;
}

static void gt_wtree_arguments_delete(void *tool_arguments)
{
  GtToolbox *wtree_toolbox = tool_arguments;
  if (wtree_toolbox != NULL) {
    gt_toolbox_delete(wtree_toolbox);
  }
}

static GtOptionParser* gt_wtree_option_parser_new(void *tool_arguments)
{
  GtToolbox *wtree_toolbox = tool_arguments;
  GtOptionParser *op;
  gt_assert(wtree_toolbox);
  op = gt_option_parser_new("[option ...] tool [argument ...]",
                            "Call an wtree manipulation tool "
                            "and pass argument(s) to it.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, wtree_toolbox);
  gt_option_parser_set_min_args(op, 1U);
  return op;
}

static int gt_wtree_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtToolbox *wtree_toolbox = tool_arguments;
  GtTool *tool = NULL;
  int had_err = 0;
  char **nargv = NULL;

  gt_error_check(err);
  gt_assert(wtree_toolbox);

  /* get wtree tools */
  tool = gt_toolbox_get_tool(wtree_toolbox, argv[parsed_args]);
  if (!tool) {
    gt_error_set(err, "wtree tool '%s' not found; option -help "
                      "lists possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call wtree tool */
  if (!had_err) {
    gt_assert(tool);
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv, err);
  }

  /* free */
  gt_cstr_array_delete(nargv);

  return had_err;
}

GtTool* gt_wtree(void)
{
  return gt_tool_new(gt_wtree_arguments_new,
                     gt_wtree_arguments_delete,
                     gt_wtree_option_parser_new,
                     NULL,
                     gt_wtree_runner);
}
