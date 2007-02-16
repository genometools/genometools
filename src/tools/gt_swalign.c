/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

#define DEFAULT_INDELSCORE -3

static int parse_options(int *indelscore, int argc, char **argv)
{
  OptionParser *op;
  Option *o;
  int parsed_args;
  op = option_parser_new("[option ...] scorematrix seq_file_1 seq_file_2",
                         "Locally align each sequence in seq_file_1 "
                         "with each sequence in seq_file_2.");
  o = option_new_int("indelscore", "set the score used for "
                     "insertions/deletions", indelscore, DEFAULT_INDELSCORE);
  option_parser_add_option(op, o);
  option_parser_parse_min_max_args(op, &parsed_args, argc, argv, versionfunc, 3,
                                   3);
  option_parser_free(op);
  return parsed_args;
}

int gt_swalign(int argc, char *argv[])
{
  Bioseq *bioseq_1 = NULL, *bioseq_2 = NULL;
  ScoreFunction *scorefunction = NULL;
  ScoreMatrix *scorematrix;
  unsigned long i, j;
  int parsed_args, indelscore;
  Alignment *a;
  Error *err = error_new();

  /* option parsing */
  parsed_args = parse_options(&indelscore, argc, argv);
  assert(parsed_args+2 < argc);

  /* init */
  /* XXX: make this more flexible */
  scorematrix  = scorematrix_read_protein(argv[parsed_args], err);
  if (scorematrix) {
    scorefunction = scorefunction_new(scorematrix, indelscore, indelscore);
    bioseq_1 = bioseq_new(argv[parsed_args+1]);
    bioseq_2 = bioseq_new(argv[parsed_args+2]);

    /* aligning all sequence combinations */
    for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
      for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
        a = swalign(bioseq_get_seq(bioseq_1, i), bioseq_get_seq(bioseq_2, j),
                    scorefunction);
        if (a) {
          alignment_show(a, stdout);
          xputchar('\n');
          alignment_free(a);
        }
      }
    }
  }

  /* free */
  bioseq_free(bioseq_2);
  bioseq_free(bioseq_1);
  scorefunction_free(scorefunction);
  error_abort(err);
  error_free(err);

  return EXIT_SUCCESS;
}
