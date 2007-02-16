/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(int *parsed_args, unsigned int *q, int argc,
                            char **argv, Error *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] seq_file_1 seq_file_2",
                         "Compute q-gram distance for each sequence "
                         "combination.");
  o = option_new_uint_min("q", "set q", q, 3, 1);
  option_parser_add_option(op, o);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, err);
  option_parser_free(op);
  return oprval;
}

int gt_qgramdist(int argc, char *argv[], Error *err)
{
  Bioseq *bioseq_1, *bioseq_2;
  unsigned long i, j, dist;
  Seq *seq_1, *seq_2;
  int parsed_args;
  unsigned int q;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &q, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args+1 < argc);

  /* make sure seq_file_1 exists */
  if (!file_exists(argv[parsed_args]))
    error("seq_file_1 \"%s\" does not exist", argv[parsed_args]);

  /* make sure seq_file_2 exists */
  if (!file_exists(argv[parsed_args+1]))
    error("seq_file_2 \"%s\" does not exist", argv[parsed_args+1]);

  /* init */
  bioseq_1 = bioseq_new(argv[parsed_args]);
  bioseq_2 = bioseq_new(argv[parsed_args+1]);

  /* compute q-gram distance for all sequence combinations */
  for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
    for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
      seq_1 = bioseq_get_seq(bioseq_1, i);
      seq_2 = bioseq_get_seq(bioseq_2, j);
      dist = qgramdist(seq_1, seq_2, q);
      printf("qgramdist_%u_(", q);
      cstr_show(seq_get_orig(seq_1), seq_length(seq_1), stdout);
      xputchar(',');
      cstr_show(seq_get_orig(seq_2), seq_length(seq_2), stdout);
      printf(")=%lu\n", dist);
    }
  }

  /* free */
  bioseq_free(bioseq_2);
  bioseq_free(bioseq_1);

  return 0;
}
