/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

#define DEFAULT_INDELSCORE -3

static OPrval parse_options(int *parsed_args, int *indelscore, int argc,
                            const char **argv, Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] scorematrix seq_file_1 seq_file_2",
                         "Locally align each sequence in seq_file_1 "
                         "with each sequence in seq_file_2.", env);
  o = option_new_int("indelscore", "set the score used for "
                     "insertions/deletions", indelscore, DEFAULT_INDELSCORE,
                     env);
  option_parser_add_option(op, o, env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 3, 3, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_swalign(int argc, const char **argv, Env *env)
{
  Bioseq *bioseq_1 = NULL, *bioseq_2 = NULL;
  ScoreFunction *scorefunction = NULL;
  ScoreMatrix *scorematrix;
  unsigned long i, j;
  int parsed_args, indelscore, had_err = 0;
  Alignment *a;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &indelscore, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args+2 < argc);

  /* init */
  /* XXX: make this more flexible */
  scorematrix  = scorematrix_read_protein(argv[parsed_args], env);
  if (scorematrix) {
    scorefunction = scorefunction_new(scorematrix, indelscore, indelscore, env);
    bioseq_1 = bioseq_new(argv[parsed_args+1], env);
    if (!bioseq_1)
      had_err = -1;
    if (!had_err) {
      bioseq_2 = bioseq_new(argv[parsed_args+2], env);
      if (!bioseq_2)
        had_err = -1;
    }

    if (!had_err) {
      /* aligning all sequence combinations */
      for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
        for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
          a = swalign(bioseq_get_seq(bioseq_1, i, env),
                      bioseq_get_seq(bioseq_2, j, env), scorefunction, env);
          if (a) {
            alignment_show(a, stdout);
            xputchar('\n');
            alignment_delete(a, env);
          }
        }
      }
    }
  }

  /* free */
  bioseq_delete(bioseq_2, env);
  bioseq_delete(bioseq_1, env);
  scorefunction_delete(scorefunction, env);

  return had_err;
}
