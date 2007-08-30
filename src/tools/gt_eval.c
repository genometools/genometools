/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gtdata.h"
#include "libgtext/stream_evaluator.h"

typedef struct {
  bool verbose,
       exondiff,
       evalLTR;
  unsigned long LTRdelta;
} EvalArguments;

static OPrval parse_options(int *parsed_args, EvalArguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *option, *ltroption, *ltrdeltaoption;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("reality_file prediction_file ",
                         "Evaluate a gene prediction against a given "
                         "``reality'' file (both in GFF3).", env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* -exondiff */
  option = option_new_bool("exondiff", "show a diff for the exons",
                           &arguments->exondiff, false, env);
  option_is_development_option(option);
  option_parser_add_option(op, option, env);

  /* -ltr */
  ltroption = option_new_bool("ltr", "evaluate a LTR retrotransposon "
                              "prediction instead of a gene prediction\n"
                              "(all LTR_retrotransposon elements are "
                              "considered to have an undetermined strand)",
                              &arguments->evalLTR, false, env);
  option_parser_add_option(op, ltroption, env);

  /* -ltrdelta */
  ltrdeltaoption = option_new_ulong("ltrdelta", "set allowed delta for LTR "
                                    "borders to be considered equal",
                                    &arguments->LTRdelta, 20, env);
  option_parser_add_option(op, ltrdeltaoption, env);

  /* option implications */
  option_imply(ltrdeltaoption, ltroption, env);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_eval(int argc, const char **argv, Env *env)
{
  GenomeStream *reality_stream,
               *prediction_stream;
  StreamEvaluator *evaluator;
  EvalArguments arguments;
  int had_err, parsed_args;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create the reality stream */
  reality_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose, env);

  /* create the prediction stream */
  prediction_stream = gff3_in_stream_new_sorted(argv[parsed_args + 1],
                                                arguments.verbose, env);

  /* create the stream evaluator */
  evaluator = stream_evaluator_new(reality_stream, prediction_stream,
                                   arguments.evalLTR, arguments.LTRdelta, env);

  /* compute the evaluation */
  had_err = stream_evaluator_evaluate(evaluator, arguments.verbose,
                                      arguments.exondiff, NULL, env);

  /* show the evaluation */
  if (!had_err)
    stream_evaluator_show(evaluator, stdout);

  /* free */
  stream_evaluator_delete(evaluator, env);
  genome_stream_delete(prediction_stream, env);
  genome_stream_delete(reality_stream, env);

  return had_err;
}
