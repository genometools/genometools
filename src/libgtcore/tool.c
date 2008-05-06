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

#include "libgtcore/ma.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/tool.h"

struct Tool {
  ToolArgumentsNew tool_arguments_new;
  ToolArgumentsDelete tool_arguments_delete;
  ToolOptionParserNew tool_option_parser_new;
  ToolArgumentsCheck tool_arguments_check;
  ToolRunner tool_runner;
};

Tool* tool_new(ToolArgumentsNew tool_arguments_new,
               ToolArgumentsDelete tool_arguments_delete,
               ToolOptionParserNew tool_option_parser_new,
               ToolArgumentsCheck tool_arguments_check,
               ToolRunner tool_runner)
{
  Tool *tool;
  assert(tool_option_parser_new && tool_runner); /* required arguments */
  /* <tool_arguments_new> and <tool_arguments_delete> imply each other */
  assert(( tool_arguments_new &&  tool_arguments_delete) ||
         (!tool_arguments_new && !tool_arguments_delete));
  tool = ma_malloc(sizeof *tool);
  tool->tool_arguments_new = tool_arguments_new;
  tool->tool_arguments_delete = tool_arguments_delete;
  tool->tool_option_parser_new = tool_option_parser_new;
  tool->tool_arguments_check = tool_arguments_check;
  tool->tool_runner = tool_runner;
  return tool;
}

int tool_run(Tool *tool, int argc, const char **argv, Error *err)
{
  void *tool_arguments = NULL;
  OptionParser *op;
  OPrval oprval;
  int parsed_args, had_err = 0;
  error_check(err);
  assert(tool);

  /* create tool arguments object */
  if (tool->tool_arguments_new)
    tool_arguments = tool->tool_arguments_new();

  /* parse options */
  op = tool->tool_option_parser_new(tool_arguments);
  oprval = option_parser_parse(op, &parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  switch (oprval) {
    case OPTIONPARSER_OK:
      break;
    case OPTIONPARSER_ERROR:
      had_err = -1;
      break;
    case OPTIONPARSER_REQUESTS_EXIT:
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

void  tool_delete(Tool *tool)
{
  if (!tool) return;
  ma_free(tool);
}
