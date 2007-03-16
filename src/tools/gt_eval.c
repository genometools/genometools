/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool verbose,
       exondiff;
} EvalArguments;

static OPrval parse_options(int *parsed_args, EvalArguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
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

  /* parse */
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
  int has_err, parsed_args;
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
  evaluator = stream_evaluator_new(reality_stream, prediction_stream, env);

  /* compute the evaluation */
  has_err = stream_evaluator_evaluate(evaluator, arguments.verbose,
                                      arguments.exondiff, env);

  /* show the evaluation */
  if (!has_err)
    stream_evaluator_show(evaluator, stdout);

  /* free */
  stream_evaluator_delete(evaluator, env);

  return has_err;
}
