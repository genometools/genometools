/*
  Copyright (c) 2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/dot_visitor.h"
#include "extended/gff3_in_stream.h"
#include "extended/visitor_stream.h"
#include "tools/gt_dot.h"

typedef struct {
} GtDotArguments;

static void* gt_dot_arguments_new(void)
{
  GtDotArguments *arguments = gt_calloc(1, sizeof *arguments);
  return arguments;
}

static void gt_dot_arguments_delete(void *tool_arguments)
{
  GtDotArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static GtOptionParser* gt_dot_option_parser_new(void *tool_arguments)
{
  GT_UNUSED GtDotArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GT_UNUSED GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [files]",
                            "Prints feature graphs in dotfile format.");

  return op;
}

static int gt_dot_runner(int argc, const char **argv, int parsed_args,
                         void *tool_arguments, GtError *err)
{
  GT_UNUSED GtDotArguments *arguments = tool_arguments;
  GtNodeStream *vs = NULL, *ins = NULL;
  GtNodeVisitor *dv = NULL;
  int had_err = 0;
  gt_error_check(err);

  ins = gt_gff3_in_stream_new_unsorted(argc - parsed_args, argv + parsed_args);
  if (ins) {
    vs = gt_visitor_stream_new(ins, (dv = gt_dot_visitor_new()));
  }
  if (vs) {
    had_err = gt_node_stream_pull(vs, err);
    gt_dot_visitor_finalize(dv);
  }
  gt_node_stream_delete(ins);
  gt_node_stream_delete(vs);
  return had_err;
}

GtTool* gt_dot(void)
{
  return gt_tool_new(gt_dot_arguments_new,
                     gt_dot_arguments_delete,
                     gt_dot_option_parser_new,
                     NULL,
                     gt_dot_runner);
}
