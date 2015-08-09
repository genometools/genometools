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
#include "tools/gt_condenseq.h"
#include "tools/gt_condenseq_compress.h"
#include "tools/gt_condenseq_extract.h"
#include "tools/gt_condenseq_info.h"
#include "tools/gt_condenseq_search.h"

static void* gt_condenseq_toolbox_new(void)
{
  GtToolbox *condenseq_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(condenseq_toolbox, "compress", gt_condenseq_compress());
  gt_toolbox_add_tool(condenseq_toolbox, "search", gt_condenseq_search());
  gt_toolbox_add_tool(condenseq_toolbox, "extract", gt_condenseq_extract());
  gt_toolbox_add_tool(condenseq_toolbox, "info", gt_condenseq_info());
  return condenseq_toolbox;
}

static void gt_condenseq_toolbox_delete(void *toolbox)
{
  GtToolbox *condenseq_toolbox = toolbox;
  gt_toolbox_delete(condenseq_toolbox);
}

static GtOptionParser* gt_condenseq_option_parser_new(void *tool_arguments)
{
  GtToolbox *condenseq_toolbox = tool_arguments;
  GtOptionParser *op;
  gt_assert(condenseq_toolbox);

  /* init */
  op = gt_option_parser_new("tool [option ...]",
                            "Call one of the CONDENSER tools to prepare "
                            "or manipulate redundancy compressed genomic data."
                            );
  gt_option_parser_set_comment_func(op, gt_toolbox_show, condenseq_toolbox);
  gt_option_parser_set_min_args(op, 1U);
  return op;
}

static int gt_condenseq_runner(int argc, const char **argv, int parsed_args,
                               void *toolbox, GtError *err)
{
  GtToolbox *condenseq_toolbox = toolbox;
  GtTool *tool;
  int had_err = 0;
  char **nargv = NULL;

  gt_error_check(err);
  gt_assert(condenseq_toolbox);

  /* get condenseq tools */
  tool = gt_toolbox_get_tool(condenseq_toolbox, argv[parsed_args]);
  if (!tool) {
    gt_error_set(err, "condenseq tool '%s' not found; option -help lists "
                      "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call condenseq tool */
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

GtTool* gt_condenseq(void)
{
  return gt_tool_new(gt_condenseq_toolbox_new,
                     gt_condenseq_toolbox_delete,
                     gt_condenseq_option_parser_new,
                     NULL,
                     gt_condenseq_runner);
}
