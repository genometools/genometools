/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "libgtcore/bioseq.h"
#include "libgtcore/fasta.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/mutate.h"
#include "tools/gt_mutate.h"

typedef struct {
  unsigned int rate; /* the mutate rate */
} MutateArguments;

static OPrval parse_options(int *parsed_args, MutateArguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] sequence_file [...]",
                         "Mutate the sequences of the given sequence_file(s) "
                         "and show them on stdout.");
  /* -rate */
  o = option_new_uint_max("rate", "set the mutation rate", &arguments->rate, 1,
                          100);
  option_parser_add_option(op, o);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  option_parser_set_min_args(op, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_mutate(int argc, const char **argv, Error *err)
{
  MutateArguments arguments;
  Bioseq *bioseq;
  unsigned long i;
  Seq *mutated_seq;
  int parsed_args, had_err = 0;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args < argc);

  while (!had_err && parsed_args < argc) {
    bioseq = bioseq_new(argv[parsed_args], err);
    if (!bioseq)
      had_err = -1;
    if (!had_err) {
      for (i = 0; i < bioseq_number_of_sequences(bioseq); i++) {
        mutated_seq = mutate(bioseq_get_description(bioseq, i),
                             bioseq_get_sequence(bioseq, i),
                             bioseq_get_sequence_length(bioseq, i),
                             bioseq_get_alpha(bioseq), arguments.rate);
        fasta_show_entry(seq_get_description(mutated_seq),
                         seq_get_orig(mutated_seq), seq_length(mutated_seq), 0);
        seq_delete(mutated_seq);
      }
    }
    bioseq_delete(bioseq);
    parsed_args++;
  }

  return had_err;
}
