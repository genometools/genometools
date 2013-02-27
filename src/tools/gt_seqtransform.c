/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/trans_table.h"
#include "tools/gt_seqtransform.h"

typedef struct {
  bool addstopaminos;
  unsigned long width;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} SeqtransformArguments;

static void* gt_seqtransform_arguments_new(void)
{
  SeqtransformArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_seqtransform_arguments_delete(void *tool_arguments)
{
  SeqtransformArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_seqtransform_option_parser_new(void *tool_arguments)
{
  SeqtransformArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o;
  gt_assert(arguments);
  op = gt_option_parser_new("[option ...] [sequence_file ...]",
                            "Perform simple transformations on the given "
                            "sequence file(s).");

  /* -addstopaminos */
  o = gt_option_new_bool("addstopaminos", "append stop amino acids ('"
                         GT_STOP_AMINO_CSTR"') to given protein sequences, if "
                         "not already present", &arguments->addstopaminos,
                         false);
  gt_option_parser_add_option(op, o);

  /* -width */
  o = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, o);

  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  return op;
}

static int gt_seqtransform_runner(int argc, const char **argv, int parsed_args,
                            void *tool_arguments, GtError *err)
{
  SeqtransformArguments *arguments = tool_arguments;
  GtBioseqIterator *bsi;
  unsigned long i;
  GtBioseq *bioseq;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  bsi = gt_bioseq_iterator_new(argc - parsed_args, argv + parsed_args);

  while (!(had_err = gt_bioseq_iterator_next(bsi, &bioseq, err)) && bioseq) {
    GtAlphabet *alphabet;
    bool is_protein;
    alphabet = gt_bioseq_get_alphabet(bioseq);
    is_protein = gt_alphabet_is_protein(alphabet);
    for (i = 0; i < gt_bioseq_number_of_sequences(bioseq); i++) {
      const char *desc, *suffix = NULL;
      char *seq;
      unsigned long seqlen;
      desc = gt_bioseq_get_description(bioseq, i);
      seq = gt_bioseq_get_sequence(bioseq, i);
      seqlen = gt_bioseq_get_sequence_length(bioseq, i);
      if (arguments->addstopaminos && is_protein && seqlen &&
          seq[seqlen-1] != GT_STOP_AMINO) {
        suffix = GT_STOP_AMINO_CSTR;
      }
      gt_fasta_show_entry_with_suffix(desc, seq, seqlen, suffix,
                                      arguments->width, arguments->outfp);
      gt_free(seq);
    }
    gt_bioseq_delete(bioseq);
  }

  gt_bioseq_iterator_delete(bsi);

  return had_err;
}

GtTool* gt_seqtransform(void)
{
  return gt_tool_new(gt_seqtransform_arguments_new,
                     gt_seqtransform_arguments_delete,
                     gt_seqtransform_option_parser_new,
                     NULL,
                     gt_seqtransform_runner);
}
