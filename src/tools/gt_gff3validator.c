/*
  Copyright (c) 2008-2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/typecheck_info.h"
#include "tools/gt_gff3validator.h"

typedef struct {
  GtTypecheckInfo *tci;
  bool strict;
} GFF3ValidatorArguments;

static void* gt_gff3validator_arguments_new(void)
{
  GFF3ValidatorArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->tci = gt_typecheck_info_new();
  return arguments;
}

static void gt_gff3validator_arguments_delete(void *tool_arguments)
{
  GFF3ValidatorArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_typecheck_info_delete(arguments->tci);
  gt_free(arguments);
}

static GtOptionParser* gt_gff3validator_option_parser_new(void *tool_arguments)
{
  GFF3ValidatorArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                         "Strictly validate given GFF3 files.");

  /* typecheck options */
  gt_typecheck_info_register_options(arguments->tci, op);

  /* -strict */
  option = gt_option_new_bool("strict", "be very strict during GFF3 parsing "
                              "(stricter than the specification requires)",
                              &arguments->strict, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_gff3validator_runner(int argc, const char **argv, int parsed_args,
                                   void *tool_arguments, GtError *err)
{
  GFF3ValidatorArguments *arguments = tool_arguments;
  GtTypeChecker *type_checker = NULL;
  GtNodeStream *gff3_in_stream;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a GFF3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream*) gff3_in_stream);

  /* set different type checker if necessary */
  if (gt_typecheck_info_option_used(arguments->tci)) {
    type_checker = gt_typecheck_info_create_type_checker(arguments->tci, err);
    if (!type_checker)
      had_err = -1;
    if (!had_err)
      gt_gff3_in_stream_set_type_checker(gff3_in_stream, type_checker);
  }

  /* enable strict mode (if necessary) */
  if (!had_err && arguments->strict)
    gt_gff3_in_stream_enable_strict_mode((GtGFF3InStream*) gff3_in_stream);

  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(gff3_in_stream, err);

  if (!had_err)
    printf("input is valid GFF3\n");

  /* free */
  gt_node_stream_delete(gff3_in_stream);
  gt_type_checker_delete(type_checker);

  return had_err;
}

GtTool* gt_gff3validator(void)
{
  return gt_tool_new(gt_gff3validator_arguments_new,
                  gt_gff3validator_arguments_delete,
                  gt_gff3validator_option_parser_new,
                  NULL,
                  gt_gff3validator_runner);
}
