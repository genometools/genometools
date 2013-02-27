/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "extended/gff3_in_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/stream_evaluator.h"
#include "tools/gt_eval.h"

typedef struct {
  bool verbose,
       exondiff,
       exondiffcollapsed,
       nuceval,
       evalLTR;
  unsigned long LTRdelta;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} EvalArguments;

static void* gt_eval_arguments_new(void)
{
  EvalArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_eval_arguments_delete(void *tool_arguments)
{
  EvalArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_eval_option_parser_new(void *tool_arguments)
{
  EvalArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *ltroption, *ltrdeltaoption;
  gt_assert(arguments);

  op = gt_option_parser_new("reference_file prediction_file ",
                            "Compare annotation files and show "
                            "accuracy measures (prediction vs. reference).");

  /* -exondiff */
  option = gt_option_new_bool("exondiff", "show a diff for the exons",
                              &arguments->exondiff, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -exondiffcollapsed */
  option = gt_option_new_bool("exondiffcollapsed", "show a diff for the "
                              "collapsed exons", &arguments->exondiffcollapsed,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -nuc */
  option = gt_option_new_bool("nuc",
                              "evaluate nucleotide level (memory consumption "
                              "is proportional to the input file sizes)",
                              &arguments->nuceval, true);
  gt_option_parser_add_option(op, option);

  /* -ltr */
  ltroption = gt_option_new_bool("ltr", "evaluate a LTR retrotransposon "
                                 "prediction instead of a gene prediction\n"
                                 "(all LTR_retrotransposon elements are "
                                 "considered to have an undetermined strand)",
                                 &arguments->evalLTR, false);
  gt_option_parser_add_option(op, ltroption);

  /* -ltrdelta */
  ltrdeltaoption = gt_option_new_ulong("ltrdelta", "set allowed delta for LTR "
                                       "borders to be considered equal",
                                       &arguments->LTRdelta, 20);
  gt_option_parser_add_option(op, ltrdeltaoption);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  /* option implications */
  gt_option_imply(ltrdeltaoption, ltroption);

  /* set comment function */
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);

  /* set minimum and maximum number of arguments */
  gt_option_parser_set_min_max_args(op, 2, 2);

  return op;
}

int gt_eval_runner(GT_UNUSED int argc, const char **argv, int parsed_args,
                   void *tool_arguments, GtError *err)
{
  EvalArguments *arguments = tool_arguments;
  GtNodeStream *reference_stream,
               *prediction_stream;
  GtStreamEvaluator *evaluator;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  /* create the reference stream */
  reference_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);
  if (arguments->verbose)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) reference_stream);

  /* create the prediction stream */
  prediction_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args + 1]);
  if (arguments->verbose)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) prediction_stream);

  /* create the stream evaluator */
  evaluator = gt_stream_evaluator_new(reference_stream, prediction_stream,
                                      arguments->nuceval, arguments->evalLTR,
                                      arguments->LTRdelta);

  /* compute the evaluation */
  had_err = gt_stream_evaluator_evaluate(evaluator, arguments->verbose,
                                         arguments->exondiff,
                                         arguments->exondiffcollapsed, NULL,
                                         err);

  /* show the evaluation */
  if (!had_err)
    gt_stream_evaluator_show(evaluator, arguments->outfp);

  /* free */
  gt_stream_evaluator_delete(evaluator);
  gt_node_stream_delete(prediction_stream);
  gt_node_stream_delete(reference_stream);

  return had_err;
}

GtTool* gt_eval(void)
{
  return gt_tool_new(gt_eval_arguments_new,
                     gt_eval_arguments_delete,
                     gt_eval_option_parser_new,
                     NULL,
                     gt_eval_runner);
}
