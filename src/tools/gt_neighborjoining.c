/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("sequence_file|example",
                         "Compute and show Neighbor-Joining tree for the "
                         "sequences in sequence file (using\nthe unit cost "
                         "edit distance as distance function). If 'example' is "
                         "given as\nsequence_file, a builtin example is used.",
                         env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

static double distfunc(unsigned long i, unsigned long j, void *data, Env *env)
{
  Bioseq *bioseq= (Bioseq*) data;
  return linearedist(bioseq_get_sequence(bioseq, i),
                     bioseq_get_sequence_length(bioseq, i),
                     bioseq_get_sequence(bioseq, j),
                     bioseq_get_sequence_length(bioseq, j), env);
}

static double exampledistfunc(unsigned long i, unsigned long j,
                              /*@unused@*/ void *data,
                              /*@unused@*/ Env *env)
{
  static const double exampledistances[5][5] =
    { {0.0   , 0.1715, 0.2147, 0.3091, 0.2326},
      {0.1715, 0.0   , 0.2991, 0.3399, 0.2058},
      {0.2147, 0.2991, 0.0   , 0.2795, 0.3943},
      {0.3091, 0.3399, 0.2795, 0.0   , 0.4289},
      {0.2326, 0.2058, 0.3943, 0.4289, 0.0   } };
  return exampledistances[i][j];
}

int gt_neighborjoining(int argc, const char **argv, Env *env)
{
  bool use_hard_coded_example = false;
  Bioseq *bioseq = NULL;
  NeighborJoining *nj = NULL;
  int parsed_args, has_err = 0;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  if (!strcmp(argv[1], "example"))
    use_hard_coded_example = true;

  if (use_hard_coded_example)
    nj = neighborjoining_new(5, NULL, exampledistfunc, env);
  else {
    bioseq = bioseq_new(argv[1], env);
    if (!bioseq)
      has_err = -1;
    if (!has_err) {
      nj = neighborjoining_new(bioseq_number_of_sequences(bioseq), bioseq,
                               distfunc, env);
    }
  }

  if (!has_err)
    neighborjoining_show_tree(nj, stdout);

  bioseq_delete(bioseq, env);
  neighborjoining_delete(nj, env);

  return has_err;
}
