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

static int parse_options(EvalArguments *arguments, int argc, char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *option;
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

  /* parse */
  parsed_args = option_parser_parse_min_max_args(op, argc, argv, versionfunc, 2,
                                                 2);
  option_parser_free(op);

  return parsed_args;
}

int gt_eval(int argc, char *argv[])
{
  Genome_stream *reality_stream,
                *prediction_stream;
  Stream_evaluator *evaluator;
  EvalArguments arguments;
  int parsed_args;

  /* option parsing */
  parsed_args = parse_options(&arguments, argc, argv);

  /* create the reality stream */
  reality_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose);

  /* create the prediction stream */
  prediction_stream = gff3_in_stream_new_sorted(argv[parsed_args + 1],
                                                arguments.verbose);

  /* create the stream evaluator */
  evaluator = stream_evaluator_new(reality_stream, prediction_stream);

  /* compute the evaluation */
  stream_evaluator_evaluate(evaluator, arguments.verbose, arguments.exondiff);

  /* show the evaluation */
  stream_evaluator_show(evaluator, stdout);

  /* free */
  stream_evaluator_free(evaluator);

  return EXIT_SUCCESS;
}
