/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  unsigned int rate; /* the mutate rate */
} MutateArguments;

static OPrval parse_options(int *parsed_args, MutateArguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] sequence_file [...]",
                         "Mutate the sequences of the given sequence_file(s) "
                         "and show them on stdout.", env);
  /* -rate */
  o = option_new_uint_max("rate", "set the mutation rate", &arguments->rate, 1,
                          100, env);
  option_parser_add_option(op, o, env);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_mutate(int argc, const char **argv, Env *env)
{
  MutateArguments arguments;
  Bioseq *bioseq;
  unsigned long i;
  Seq *mutated_seq;
  int parsed_args, has_err = 0;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args < argc);

  while (!has_err && parsed_args < argc) {
    bioseq = bioseq_new(argv[parsed_args], env);
    if (!bioseq)
      has_err = -1;
    if (!has_err) {
      for (i = 0; i < bioseq_number_of_sequences(bioseq); i++) {
        mutated_seq = mutate(bioseq_get_description(bioseq, i),
                             bioseq_get_sequence(bioseq, i),
                             bioseq_get_sequence_length(bioseq, i),
                             bioseq_get_alpha(bioseq, env), arguments.rate,
                             env);
        fasta_show_entry(seq_get_description(mutated_seq),
                         seq_get_orig(mutated_seq), seq_length(mutated_seq), 0);
        seq_delete(mutated_seq, env);
      }
    }
    bioseq_delete(bioseq, env);
    parsed_args++;
  }

  return has_err;
}
