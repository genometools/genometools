/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/versionfunc.h"
#include "core/tool_api.h"

struct GtTool {
  GtToolArgumentsNew tool_arguments_new;
  GtToolArgumentsDelete tool_arguments_delete;
  GtToolOptionParserNew tool_option_parser_new;
  GtToolArgumentsCheck tool_arguments_check;
  GtToolRunner tool_runner;
};

GtTool* gt_tool_new(GtToolArgumentsNew tool_arguments_new,
                    GtToolArgumentsDelete tool_arguments_delete,
                    GtToolOptionParserNew tool_option_parser_new,
                    GtToolArgumentsCheck tool_arguments_check,
                    GtToolRunner tool_runner)
{
  GtTool *tool;
  gt_assert(tool_option_parser_new && tool_runner); /* required arguments */
  /* <tool_arguments_new> and <tool_arguments_delete> imply each other */
  gt_assert(( tool_arguments_new &&  tool_arguments_delete) ||
         (!tool_arguments_new && !tool_arguments_delete));
  tool = gt_malloc(sizeof *tool);
  tool->tool_arguments_new = tool_arguments_new;
  tool->tool_arguments_delete = tool_arguments_delete;
  tool->tool_option_parser_new = tool_option_parser_new;
  tool->tool_arguments_check = tool_arguments_check;
  tool->tool_runner = tool_runner;
  return tool;
}

int gt_tool_run(GtTool *tool, int argc, const char **argv, GtError *err)
{
  void *tool_arguments = NULL;
  GtOptionParser *op;
  GtOPrval oprval;
  int parsed_args, had_err = 0;
  gt_error_check(err);
  gt_assert(tool);

  /* create tool arguments object */
  if (tool->tool_arguments_new)
    tool_arguments = tool->tool_arguments_new();

  /* parse options */
  op = tool->tool_option_parser_new(tool_arguments);
  oprval = gt_option_parser_parse(op, &parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  switch (oprval) {
    case GT_OPTION_PARSER_OK:
      break;
    case GT_OPTION_PARSER_ERROR:
      had_err = -1;
      break;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
      if (tool_arguments)
        tool->tool_arguments_delete(tool_arguments);
      return 0;
  }

  /* check tool arguments */
  if (!had_err && tool->tool_arguments_check) {
    had_err = tool->tool_arguments_check(argc - parsed_args , tool_arguments,
                                         err);
  }

  /* run tool */
  if (!had_err)
    had_err = tool->tool_runner(argc, argv, parsed_args, tool_arguments, err);

  /* delete tool argument object */
  if (tool_arguments)
    tool->tool_arguments_delete(tool_arguments);

  /* return */
  if (had_err)
    return -1;
  return 0;
}

void  gt_tool_delete(GtTool *tool)
{
  if (!tool) return;
  gt_free(tool);
}
