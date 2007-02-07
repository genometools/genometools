/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  unsigned int show,
               consensus,
               sumofpairs;
} MSAparse_arguments;

static int parse_options(MSAparse_arguments *arguments, int argc, char **argv)
{
  OptionParser *op;
  Option *o;
  int parsed_args;
  op = option_parser_new("[option ...] MSA_file",
                         "Parse multiple sequence alignment (MSA) file and "
                         "optionally show score(s).");
  /* -show */
  o = option_new_boolean("show", "show the parsed MSA on stdout",
                         &arguments->show, 0);
  option_parser_add_option(op, o);
  /* -consensus */
  o = option_new_boolean("consensus", "show consensus distance",
                         &arguments->consensus, 0);
  option_parser_add_option(op, o);
  /* -sumofpairs */
  o = option_new_boolean("sumofpairs", "show optimal sum of pairwise scores",
                         &arguments->sumofpairs, 0);
  option_parser_add_option(op, o);
  /* parse */
  parsed_args = option_parser_parse_min_max_args(op, argc, argv, versionfunc, 1,
                                                 1);
  option_parser_free(op);
  return parsed_args;
}

int gt_msaparse(int argc, char *argv[])
{
  MSAparse_arguments arguments;
  int parsed_args;
  MSA *msa;

  /* option parsing */
  parsed_args = parse_options(&arguments, argc, argv);

  /* make sure sequence_file exists */
  assert(parsed_args < argc);
  if (!file_exists(argv[parsed_args]))
    error("MSA_file '%s' does not exist", argv[parsed_args]);

  /* multiple sequence alignment construction */
  msa = msa_new(argv[parsed_args]);

  /* output */
  if (arguments.show)
    msa_show(msa);
  if (arguments.consensus)
    printf("consensus distance: %lu\n", msa_consensus_distance(msa));
  if (arguments.sumofpairs)
    printf("sum of pairwise scores: %lu\n", msa_sum_of_pairwise_scores(msa));

  /* free */
  msa_free(msa);

  return EXIT_SUCCESS;
}
