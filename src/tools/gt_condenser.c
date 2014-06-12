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
#include "tools/gt_condenser.h"
#include "tools/gt_condenser_compress.h"
#include "tools/gt_condenser_compsearch.h"
#include "tools/gt_condenser_extract.h"
#include "tools/gt_condenser_info.h"

typedef struct {
} GtCondenserArguments;

static void* gt_condenser_arguments_new(void)
{
  GtToolbox *condenser_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(condenser_toolbox, "compress", gt_condenser_compress());
  gt_toolbox_add_tool(condenser_toolbox, "compsearch",
                      gt_condenser_compsearch());
  gt_toolbox_add_tool(condenser_toolbox, "extract", gt_condenser_extract());
  gt_toolbox_add_tool(condenser_toolbox, "info", gt_condenser_info());
  return condenser_toolbox;
}

static void gt_condenser_arguments_delete(void *tool_arguments)
{
  GtToolbox *condenser_toolbox = tool_arguments;
  if (condenser_toolbox != NULL) {
    gt_toolbox_delete(condenser_toolbox);
  }
}

static GtOptionParser* gt_condenser_option_parser_new(void *tool_arguments)
{
  GtToolbox *condenser_toolbox = tool_arguments;
  GtOptionParser *op;
  gt_assert(condenser_toolbox);

  /* init */
  op = gt_option_parser_new("[option ...] [file]",
                            "Call one of the CONDENSER tools to prepare "
                            "or manipulate redundancie compressed genomic data."
                            );
  gt_option_parser_set_comment_func(op, gt_toolbox_show, condenser_toolbox);
  gt_option_parser_set_min_args(op, 1U);
  return op;
}

static int gt_condenser_runner(int argc, const char **argv, int parsed_args,
                          void *tool_arguments, GT_UNUSED GtError *err)
{
  GtToolbox *condenser_toolbox = tool_arguments;
  GtTool *tool;
  int had_err = 0;
  char **nargv = NULL;

  gt_error_check(err);
  gt_assert(condenser_toolbox);

  /* get condenser tools */
  tool = gt_toolbox_get_tool(condenser_toolbox, argv[parsed_args]);
  if (!tool) {
    gt_error_set(err, "condenser tool '%s' not found; option -help lists "
                      "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call condenser tool */
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

GtTool* gt_condenser(void)
{
  return gt_tool_new(gt_condenser_arguments_new,
                     gt_condenser_arguments_delete,
                     gt_condenser_option_parser_new,
                     NULL,
                     gt_condenser_runner);
}
