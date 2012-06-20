/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@studium.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr_array.h"
#include "core/unused_api.h"
#include "extended/toolbox.h"
#include "tools/gt_compreads.h"
#include "tools/gt_compreads_refcompress.h"
#include "tools/gt_compreads_refdecompress.h"
#include "tools/gt_compreads_compress.h"
#include "tools/gt_compreads_decompress.h"

static void* gt_compreads_arguments_new(void)
{
  GtToolbox *compreads_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(compreads_toolbox,
      "refcompress", gt_compreads_refcompress());
  gt_toolbox_add_tool(compreads_toolbox,
      "refdecompress", gt_compreads_refdecompress());
  gt_toolbox_add_tool(compreads_toolbox,
      "compress", gt_compreads_compress());
  gt_toolbox_add_tool(compreads_toolbox,
      "decompress", gt_compreads_decompress());
  return compreads_toolbox;
}

static void gt_compreads_arguments_delete(void *tool_arguments)
{
  GtToolbox *compreads_toolbox = tool_arguments;
  if (!compreads_toolbox)
    return;
  gt_toolbox_delete(compreads_toolbox);
}

static GtOptionParser* gt_compreads_option_parser_new(void *tool_arguments)
{
  GtToolbox *compreads_toolbox = tool_arguments;
  GtOptionParser *op;
  gt_assert(compreads_toolbox);
  op = gt_option_parser_new("[option ...] tool [argument ...]",
                            "Call fastq file compression tool <tool>.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, compreads_toolbox);
  gt_option_parser_set_min_args(op, 1U);
  return op;
}

static int gt_compreads_runner(GT_UNUSED int argc,
                         const char **argv,
                         int parsed_args,
                         void *tool_arguments, GtError *err)
{
  GtToolbox *compreads_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  int had_err = 0;
  char **nargv = NULL;

  gt_error_check(err);
  gt_assert(compreads_toolbox);
  if (!gt_toolbox_has_tool(compreads_toolbox, argv[parsed_args])) {
    gt_error_set(err, "compreads tool '%s' not found; option -help lists "
                      "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  if (!had_err) {
    if (!(toolfunc = gt_toolbox_get(compreads_toolbox, argv[parsed_args]))) {
      tool = gt_toolbox_get_tool(compreads_toolbox, argv[parsed_args]);
      gt_assert(tool);
    }

    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    if (toolfunc != NULL)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv,
                            err);
  }
  gt_cstr_array_delete(nargv);

  return had_err;
}

GtTool* gt_compreads(void)
{
  return gt_tool_new(gt_compreads_arguments_new,
                  gt_compreads_arguments_delete,
                  gt_compreads_option_parser_new,
                  NULL,
                  gt_compreads_runner);
}
