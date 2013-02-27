/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/bioseq_iterator.h"
#include "core/fasta.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "extended/gtdatahelp.h"
#include "extended/mutate.h"
#include "tools/gt_seqmutate.h"

typedef struct {
  unsigned int rate; /* the mutate rate */
  unsigned long width;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} MutateArguments;

static void* gt_seqmutate_arguments_new(void)
{
  MutateArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_seqmutate_arguments_delete(void *tool_arguments)
{
  MutateArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_seqmutate_option_parser_new(void *tool_arguments)
{
  MutateArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o;
  gt_assert(arguments);
  op = gt_option_parser_new("[option ...] [sequence_file ...]",
                            "Mutate the sequences of the given sequence "
                            "file(s).");
  /* -rate */
  o = gt_option_new_uint_max("rate", "set the mutation rate", &arguments->rate,
                             1, 100);
  gt_option_parser_add_option(op, o);

  /* -width */
  o = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, o);

  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);

  return op;
}

static int gt_seqmutate_runner(int argc, const char **argv, int parsed_args,
                            void *tool_arguments, GtError *err)
{
  MutateArguments *arguments = tool_arguments;
  GtBioseqIterator *bsi;
  unsigned long i;
  GtBioseq *bioseq;
  GtSeq *mutated_seq;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  bsi = gt_bioseq_iterator_new(argc - parsed_args, argv + parsed_args);

  while (!(had_err = gt_bioseq_iterator_next(bsi, &bioseq, err)) && bioseq) {
    char *seq;
    for (i = 0; i < gt_bioseq_number_of_sequences(bioseq); i++) {
      seq = gt_bioseq_get_sequence(bioseq, i);
      mutated_seq = gt_mutate_seq(gt_bioseq_get_description(bioseq, i),
                                  seq,
                                  gt_bioseq_get_sequence_length(bioseq, i),
                                  gt_bioseq_get_alphabet(bioseq),
                                  arguments->rate);
      gt_fasta_show_entry(gt_seq_get_description(mutated_seq),
                          gt_seq_get_orig(mutated_seq),
                          gt_seq_length(mutated_seq),
                          arguments->width, arguments->outfp);
      gt_free(seq);
      gt_seq_delete(mutated_seq);
    }
    gt_bioseq_delete(bioseq);
  }

  gt_bioseq_iterator_delete(bsi);

  return had_err;
}

GtTool* gt_seqmutate(void)
{
  return gt_tool_new(gt_seqmutate_arguments_new,
                     gt_seqmutate_arguments_delete,
                     gt_seqmutate_option_parser_new,
                     NULL,
                     gt_seqmutate_runner);
}
