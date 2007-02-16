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

static int parse_options(Costs *costs, int argc, char **argv)
{
  OptionParser *op;
  Option *option;
  int parsed_args;
  op = option_parser_new("[option ...] seq_file_1 seq_file_2",
                         "Globally align each sequence in seq_file_1 with each "
                         "sequence in seq_file_2 (affine gap costs).");
  option = option_new_int("rep", "set replacement cost",
                          &costs->replacement_cost, 1);
  option_parser_add_option(op, option);
  option = option_new_int("gapopen", "set gap opening cost",
                          &costs->gap_opening_cost, 3);
  option_parser_add_option(op, option);
  option = option_new_int("gapext", "set gap extension cost",
                          &costs->gap_extension_cost, 1);
  option_parser_add_option(op, option);
  option_parser_parse_min_max_args(op, &parsed_args, argc, argv, versionfunc, 2,
                                   2);
  option_parser_free(op);
  return parsed_args;
}

int gt_affinealign(int argc, char *argv[])
{
  Bioseq *bioseq_1, *bioseq_2;
  unsigned long i, j;
  int parsed_args;
  Alignment *a;
  Costs costs;

  /* option parsing */
  parsed_args = parse_options(&costs, argc, argv);
  assert(parsed_args+1 < argc);

  /* init */
  bioseq_1 = bioseq_new(argv[parsed_args]);
  bioseq_2 = bioseq_new(argv[parsed_args+1]);

  /* aligning all sequence combinations */
  for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
    for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
      a = affinealign(bioseq_get_sequence(bioseq_1, i),
                      bioseq_get_sequence_length(bioseq_1, i),
                      bioseq_get_sequence(bioseq_2, j),
                      bioseq_get_sequence_length(bioseq_2, j),
                      costs.replacement_cost, costs.gap_opening_cost,
                      costs.gap_extension_cost);
      alignment_show(a, stdout);
      xputchar('\n');
      alignment_free(a);
    }
  }

  /* free */
  bioseq_free(bioseq_2);
  bioseq_free(bioseq_1);

  return EXIT_SUCCESS;
}
