/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/option.h"
#include "core/versionfunc.h"
#include "extended/gff3_in_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/stream_evaluator.h"
#include "tools/gt_eval.h"

typedef struct {
  bool verbose,
       exondiff,
       nuceval,
       evalLTR;
  unsigned long LTRdelta;
} EvalArguments;

static GtOPrval parse_options(int *parsed_args, EvalArguments *arguments,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *option, *ltroption, *ltrdeltaoption;
  GtOPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("reality_file prediction_file ",
                         "Evaluate a gene prediction against a given "
                         "``reality'' file (both in GFF3).");

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* -exondiff */
  option = gt_option_new_bool("exondiff", "show a diff for the exons",
                           &arguments->exondiff, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -nuc */
  option = gt_option_new_bool("nuc",
                              "evaluate nucleotide level (memory intense)",
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

  /* option implications */
  gt_option_imply(ltrdeltaoption, ltroption);

  /* parse */
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_min_max_args(op, 2, 2);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_eval(int argc, const char **argv, GtError *err)
{
  GtNodeStream *reality_stream,
               *prediction_stream;
  GtStreamEvaluator *evaluator;
  EvalArguments arguments;
  int had_err, parsed_args;
  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create the reality stream */
  reality_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);
  if (arguments.verbose)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) reality_stream);

  /* create the prediction stream */
  prediction_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args + 1]);
  if (arguments.verbose)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) prediction_stream);

  /* create the stream evaluator */
  evaluator = gt_stream_evaluator_new(reality_stream, prediction_stream,
                                   arguments.nuceval, arguments.evalLTR,
                                   arguments.LTRdelta);

  /* compute the evaluation */
  had_err = gt_stream_evaluator_evaluate(evaluator, arguments.verbose,
                                      arguments.exondiff, NULL, err);

  /* show the evaluation */
  if (!had_err)
    gt_stream_evaluator_show(evaluator, stdout);

  /* free */
  gt_stream_evaluator_delete(evaluator);
  gt_node_stream_delete(prediction_stream);
  gt_node_stream_delete(reality_stream);

  return had_err;
}
