/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/fasta_separator.h"
#include "libgtcore/unused.h"
#include "libgtexercise/simple_bioseq.h"
#include "tools/gt_fastaparser.h"

static OptionParser* gt_fastaparser_option_parser_new(UNUSED
                                                      void *tool_arguments)
{
  OptionParser *op;
  op = option_parser_new("[option ...] fasta_file",
                         "Parser fasta_file and show it on stdout.");
  option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static int gt_fastaparser_runner(UNUSED int argc, const char **argv,
                                 int parsed_args, UNUSED void *tool_arguments,
                                 UNUSED Error *err)
{
  unsigned long i;
  SimpleBioseq *simple_bioseq;

  error_check(err);

  simple_bioseq = simple_bioseq_new(argv[parsed_args]);
  for (i = 0; i < simple_bioseq_number_of_sequences(simple_bioseq); i++) {
    printf("%c%s\n%s\n", FASTA_SEPARATOR,
           simple_bioseq_get_description(simple_bioseq, i),
           simple_bioseq_get_sequence(simple_bioseq, i));
    printf("sequence #%lu length: %lu\n", i,
           simple_bioseq_get_sequence_length(simple_bioseq, i));
  }
  simple_bioseq_delete(simple_bioseq);

  return 0;
}

Tool* gt_fastaparser(void)
{
  return tool_new(NULL,
                  NULL,
                  gt_fastaparser_option_parser_new,
                  NULL,
                  gt_fastaparser_runner);
}
