/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static int parse_options(unsigned int *all, int argc, char **argv)
{
  OptionParser *op;
  Option *option;
  int parsed_args;
  op = option_parser_new("[option ...] seq_file_1 seq_file_2",
                         "Globally align each sequence in seq_file_1 with each "
                         "sequence in seq_file_2.");
  option = option_new_boolean("all", "show all optimal alignments instead of "
                              "just one", all, 0);
  option_parser_add_option(op, option);
  parsed_args = option_parser_parse_min_max_args(op, argc, argv, versionfunc, 2,
                                                 2);
  option_parser_free(op);
  return parsed_args;
}

void show_alignment(const Alignment *a, void *data)
{
  assert(a && !data);
  alignment_show(a, stdout);
  xputchar('\n');
}

void show_aligns(unsigned long aligns, void *data)
{
  assert(aligns && !data);
  printf("number of optimal alignments: %lu\n\n", aligns);
}

int gt_align(int argc, char *argv[])
{
  Bioseq *bioseq_1, *bioseq_2;
  unsigned long i, j;
  unsigned int all = 0;
  int parsed_args;
  Alignment *a;

  /* option parsing */
  parsed_args = parse_options(&all, argc, argv);
  assert(parsed_args+1 < argc);

  /* init */
  bioseq_1 = bioseq_new(argv[parsed_args]);
  bioseq_2 = bioseq_new(argv[parsed_args+1]);

  /* aligning all sequence combinations */
  for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
    for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
      if (all) {
        align_all(bioseq_get_sequence(bioseq_1, i),
                  bioseq_get_sequence_length(bioseq_1, i),
                  bioseq_get_sequence(bioseq_2, j),
                  bioseq_get_sequence_length(bioseq_2, j),
                  show_alignment, show_aligns, NULL);
      }
      else {
        a = align(bioseq_get_sequence(bioseq_1, i),
                  bioseq_get_sequence_length(bioseq_1, i),
                  bioseq_get_sequence(bioseq_2, j),
                  bioseq_get_sequence_length(bioseq_2, j));
        alignment_show(a, stdout);
        xputchar('\n');
        alignment_free(a);
      }
    }
  }

  /* free */
  bioseq_free(bioseq_2);
  bioseq_free(bioseq_1);

  return EXIT_SUCCESS;
}
