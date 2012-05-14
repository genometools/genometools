/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
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

#include "core/ma.h"
#include "core/unused_api.h"
#include "tools/gt_encseq2kmers.h"

typedef struct {
  unsigned int kmersize;
  GtStr  *encseqinput;
  bool verbose;
} GtEncseq2kmersArguments;

static void* gt_encseq2kmers_arguments_new(void)
{
  GtEncseq2kmersArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->encseqinput= gt_str_new();
  return arguments;
}

static void gt_encseq2kmers_arguments_delete(void *tool_arguments)
{
  GtEncseq2kmersArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq2kmers_option_parser_new(void *tool_arguments)
{
  GtEncseq2kmersArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-ii <indexname> -k <kmersize> [option ...]",
                         "Collect and process kmers of an encseq.");

  /* -k */
  option = gt_option_new_uint_min_max("k", "specify the kmer size",
                                  &arguments->kmersize, 31U, 1U, 32U);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -ii */
  option = gt_option_new_string("ii", "specify the input sequence",
                                arguments->encseqinput, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_max_args(op, 0);

  return op;
}

static int gt_encseq2kmers_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GT_UNUSED GtEncseq2kmersArguments *arguments = tool_arguments;
  int had_err = 0;
  return had_err;
}

static int gt_encseq2kmers_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    GT_UNUSED void *tool_arguments, GT_UNUSED GtError *err)
{
  GtEncseq2kmersArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

GtTool* gt_encseq2kmers(void)
{
  return gt_tool_new(gt_encseq2kmers_arguments_new,
                  gt_encseq2kmers_arguments_delete,
                  gt_encseq2kmers_option_parser_new,
                  gt_encseq2kmers_arguments_check,
                  gt_encseq2kmers_runner);
}
