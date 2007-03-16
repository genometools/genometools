/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  int replacement_cost,
      gap_opening_cost,
      gap_extension_cost;
} Costs;

static OPrval parse_options(int *parsed_args, Costs *costs, int argc,
                            const char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] seq_file_1 seq_file_2",
                         "Globally align each sequence in seq_file_1 with each "
                         "sequence in seq_file_2 (affine gap costs).", env);
  option = option_new_int("rep", "set replacement cost",
                          &costs->replacement_cost, 1, env);
  option_parser_add_option(op, option, env);
  option = option_new_int("gapopen", "set gap opening cost",
                          &costs->gap_opening_cost, 3, env);
  option_parser_add_option(op, option, env);
  option = option_new_int("gapext", "set gap extension cost",
                          &costs->gap_extension_cost, 1, env);
  option_parser_add_option(op, option, env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_affinealign(int argc, const char **argv, Env *env)
{
  Bioseq *bioseq_1, *bioseq_2 = NULL;
  unsigned long i, j;
  int parsed_args, has_err = 0;
  Alignment *a;
  Costs costs;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &costs, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args+1 < argc);

  /* init */
  bioseq_1 = bioseq_new(argv[parsed_args], env);
  if (!bioseq_1)
     has_err = -1;
  if (!has_err) {
    bioseq_2 = bioseq_new(argv[parsed_args+1], env);
    if (!bioseq_2)
      has_err = -1;
  }

  /* aligning all sequence combinations */
  if (!has_err) {
    for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
      for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
        a = affinealign(bioseq_get_sequence(bioseq_1, i),
                        bioseq_get_sequence_length(bioseq_1, i),
                        bioseq_get_sequence(bioseq_2, j),
                        bioseq_get_sequence_length(bioseq_2, j),
                        costs.replacement_cost, costs.gap_opening_cost,
                        costs.gap_extension_cost, env);
        alignment_show(a, stdout);
        xputchar('\n');
        alignment_delete(a, env);
      }
    }
  }

  /* free */
  bioseq_delete(bioseq_2, env);
  bioseq_delete(bioseq_1, env);

  return has_err;
}
