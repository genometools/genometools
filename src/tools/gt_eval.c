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

#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/stream_evaluator.h"
#include "tools/gt_eval.h"

typedef struct {
  bool verbose,
       exondiff,
       nuceval,
       evalLTR;
  unsigned long LTRdelta;
} EvalArguments;

static OPrval parse_options(int *parsed_args, EvalArguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *option, *ltroption, *ltrdeltaoption;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("reality_file prediction_file ",
                         "Evaluate a gene prediction against a given "
                         "``reality'' file (both in GFF3).");

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* -exondiff */
  option = option_new_bool("exondiff", "show a diff for the exons",
                           &arguments->exondiff, false);
  option_is_development_option(option);
  option_parser_add_option(op, option);

  /* -nuc */
  option = option_new_bool("nuc", "evaluate nucleotide level (memory intense)",
                           &arguments->nuceval, true);
  option_parser_add_option(op, option);

  /* -ltr */
  ltroption = option_new_bool("ltr", "evaluate a LTR retrotransposon "
                              "prediction instead of a gene prediction\n"
                              "(all LTR_retrotransposon elements are "
                              "considered to have an undetermined strand)",
                              &arguments->evalLTR, false);
  option_parser_add_option(op, ltroption);

  /* -ltrdelta */
  ltrdeltaoption = option_new_ulong("ltrdelta", "set allowed delta for LTR "
                                    "borders to be considered equal",
                                    &arguments->LTRdelta, 20);
  option_parser_add_option(op, ltrdeltaoption);

  /* option implications */
  option_imply(ltrdeltaoption, ltroption);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  option_parser_set_min_max_args(op, 2, 2);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_eval(int argc, const char **argv, Error *err)
{
  GenomeStream *reality_stream,
               *prediction_stream;
  StreamEvaluator *evaluator;
  EvalArguments arguments;
  int had_err, parsed_args;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create the reality stream */
  reality_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose);

  /* create the prediction stream */
  prediction_stream = gff3_in_stream_new_sorted(argv[parsed_args + 1],
                                                arguments.verbose);

  /* create the stream evaluator */
  evaluator = stream_evaluator_new(reality_stream, prediction_stream,
                                   arguments.nuceval, arguments.evalLTR,
                                   arguments.LTRdelta);

  /* compute the evaluation */
  had_err = stream_evaluator_evaluate(evaluator, arguments.verbose,
                                      arguments.exondiff, NULL, err);

  /* show the evaluation */
  if (!had_err)
    stream_evaluator_show(evaluator, stdout);

  /* free */
  stream_evaluator_delete(evaluator);
  genome_stream_delete(prediction_stream);
  genome_stream_delete(reality_stream);

  return had_err;
}
