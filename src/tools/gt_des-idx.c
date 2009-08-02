/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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
#include "core/option.h"
#include "core/unused_api.h"
#include "core/tool.h"
#include "match/giextract.h"

typedef struct {
  bool verbose;
  GtStr *indexname;
} DesIdxArguments;

static void* gt_des_idx_arguments_new(void)
{
  DesIdxArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->verbose = false;
  arguments->indexname = gt_str_new();
  return arguments;
}

static void gt_des_idx_arguments_delete(void *tool_arguments)
{
  DesIdxArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->indexname);
  gt_free(arguments);
}

static GtOptionParser* gt_des_idx_option_parser_new(void *tool_arguments)
{
  DesIdxArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] indexname",
                            "Index the keys of the form |key| in des file.");

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 1U, 1U);

  return op;
}

static int gt_des_idx_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
                             GT_UNUSED int parsed_args, void *tool_arguments,
                             GtError *err)
{
  DesIdxArguments *arguments = tool_arguments;

  gt_error_check(err);
  gt_assert(arguments);

  return gt_extractkeysfromdesfile(arguments->indexname, err);
}

GtTool* gt_des_idx(void)
{
  return gt_tool_new(gt_des_idx_arguments_new,
                     gt_des_idx_arguments_delete,
                     gt_des_idx_option_parser_new,
                     NULL,
                     gt_des_idx_runner);
}
