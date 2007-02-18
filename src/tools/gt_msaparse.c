/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool show,
       consensus,
       sumofpairs;
} MSAparse_arguments;

static OPrval parse_options(int *parsed_args, MSAparse_arguments *arguments,
                            int argc, char **argv, Error *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] MSA_file",
                         "Parse multiple sequence alignment (MSA) file and "
                         "optionally show score(s).");
  /* -show */
  o = option_new_bool("show", "show the parsed MSA on stdout", &arguments->show,
                      false);
  option_parser_add_option(op, o);
  /* -consensus */
  o = option_new_bool("consensus", "show consensus distance",
                      &arguments->consensus, false);
  option_parser_add_option(op, o);
  /* -sumofpairs */
  o = option_new_bool("sumofpairs", "show optimal sum of pairwise scores",
                      &arguments->sumofpairs, false);
  option_parser_add_option(op, o);
  /* parse */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, err);
  option_parser_free(op);
  return oprval;
}

int gt_msaparse(int argc, char *argv[], Error *err)
{
  MSAparse_arguments arguments;
  int parsed_args, has_err = 0;
  MSA *msa = NULL;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* make sure sequence_file exists */
  assert(parsed_args < argc);
  if (!file_exists(argv[parsed_args])) {
    error_set(err, "MSA_file '%s' does not exist", argv[parsed_args]);
    has_err = -1;
  }

  if (!has_err) {
    /* multiple sequence alignment construction */
    msa = msa_new(argv[parsed_args]);

    /* output */
    if (arguments.show)
      msa_show(msa);
    if (arguments.consensus)
      printf("consensus distance: %lu\n", msa_consensus_distance(msa));
    if (arguments.sumofpairs)
      printf("sum of pairwise scores: %lu\n", msa_sum_of_pairwise_scores(msa));
  }

  /* free */
  msa_free(msa);

  return has_err;
}
