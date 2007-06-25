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
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] MSA_file",
                         "Parse multiple sequence alignment (MSA) file and "
                         "optionally show score(s).", env);
  /* -show */
  o = option_new_bool("show", "show the parsed MSA on stdout", &arguments->show,
                      false, env);
  option_parser_add_option(op, o, env);
  /* -consensus */
  o = option_new_bool("consensus", "show consensus distance",
                      &arguments->consensus, false, env);
  option_parser_add_option(op, o, env);
  /* -sumofpairs */
  o = option_new_bool("sumofpairs", "show optimal sum of pairwise scores",
                      &arguments->sumofpairs, false, env);
  option_parser_add_option(op, o, env);
  /* parse */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_msaparse(int argc, const char **argv, Env *env)
{
  MSAparse_arguments arguments;
  int parsed_args, had_err = 0;
  MSA *msa = NULL;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* make sure sequence_file exists */
  assert(parsed_args < argc);
  if (!file_exists(argv[parsed_args])) {
    env_error_set(env, "MSA_file '%s' does not exist", argv[parsed_args]);
    had_err = -1;
  }

  if (!had_err) {
    /* multiple sequence alignment construction */
    msa = msa_new(argv[parsed_args], env);
    if (!msa)
      had_err = -1;
  }

  if (!had_err) {
    /* output */
    if (arguments.show)
      msa_show(msa);
    if (arguments.consensus)
      printf("consensus distance: %lu\n", msa_consensus_distance(msa, env));
    if (arguments.sumofpairs) {
      printf("sum of pairwise scores: %lu\n",
             msa_sum_of_pairwise_scores(msa, env));
    }
  }

  /* free */
  msa_delete(msa, env);

  return had_err;
}
