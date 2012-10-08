/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include "core/option_api.h"
#include "core/unused_api.h"
#include "extended/toolbox.h"
#include "extended/gtdatahelp.h"
#include "match/rdj-version.h"
#include "tools/gt_readjoiner_prefilter.h"
#include "tools/gt_readjoiner_overlap.h"
#include "tools/gt_readjoiner_assembly.h"
#include "tools/gt_readjoiner_asqg.h"
#include "tools/gt_readjoiner_cnttest.h"
#include "tools/gt_readjoiner_spmtest.h"
#include "tools/gt_readjoiner_correct.h"

static void* gt_readjoiner_arguments_new(void)
{
  GtToolbox *readjoiner_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(readjoiner_toolbox, "prefilter",
      gt_readjoiner_prefilter());
  gt_toolbox_add_tool(readjoiner_toolbox, "overlap", gt_readjoiner_overlap());
  gt_toolbox_add_tool(readjoiner_toolbox, "assembly", gt_readjoiner_assembly());
  gt_toolbox_add_hidden_tool(readjoiner_toolbox, "asqg", gt_readjoiner_asqg());
  gt_toolbox_add_hidden_tool(readjoiner_toolbox, "cnttest",
      gt_readjoiner_cnttest());
  gt_toolbox_add_hidden_tool(readjoiner_toolbox, "spmtest",
      gt_readjoiner_spmtest());
  gt_toolbox_add_hidden_tool(readjoiner_toolbox, "correct",
      gt_readjoiner_correct());
  return readjoiner_toolbox;
}

static void gt_readjoiner_arguments_delete(void *tool_arguments)
{
  GtToolbox *readjoiner_toolbox = tool_arguments;
  if (!readjoiner_toolbox) return;
  gt_toolbox_delete(readjoiner_toolbox);
}

static GtOptionParser* gt_readjoiner_option_parser_new(GT_UNUSED
    void *tool_arguments)
{
  GtOptionParser *op;
  op = gt_option_parser_new("[option ...] tool [argument ...]",
                 "Readjoiner: a string graph-based sequence assembler.");
  gt_option_parser_set_version_func(op, gt_readjoiner_show_version);
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_min_args(op, 1);
  return op;
}

static int gt_readjoiner_runner(int argc, const char **argv, int parsed_args,
                         void *tool_arguments, GtError *err)
{
  GtToolbox *readjoiner_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  int had_err = 0;
  char **nargv = NULL;

  gt_error_check(err);
  gt_assert(readjoiner_toolbox);

  /* get readjoiner tools */
  if (!gt_toolbox_has_tool(readjoiner_toolbox, argv[parsed_args])) {
    gt_error_set(err, "readjoiner tool '%s' not found; option -help lists "
                      "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call readjoiner tool */
  if (!had_err) {
    if (!(toolfunc = gt_toolbox_get(readjoiner_toolbox, argv[parsed_args]))) {
      tool = gt_toolbox_get_tool(readjoiner_toolbox, argv[parsed_args]);
      gt_assert(tool);
    }
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    if (toolfunc)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv,
                            err);
  }

  /* free */
  gt_cstr_array_delete(nargv);

  return had_err;
}

GtTool* gt_readjoiner(void)
{
  return gt_tool_new(gt_readjoiner_arguments_new,
                     gt_readjoiner_arguments_delete,
                     gt_readjoiner_option_parser_new,
                     NULL,
                     gt_readjoiner_runner);
}
