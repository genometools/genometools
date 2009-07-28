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

#include <string.h>
#include "core/bioseq_iterator.h"
#include "core/fasta.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/undef.h"
#include "tools/gt_seqfilter.h"

typedef struct {
  unsigned long minlength,
                maxlength;
} GtSeqFilterArguments;

static void* gt_seqfilter_arguments_new(void)
{
  GtSeqFilterArguments *arguments = gt_calloc(1, sizeof *arguments);
  return arguments;
}

static void gt_seqfilter_arguments_delete(void *tool_arguments)
{
  GtSeqFilterArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static GtOptionParser* gt_seqfilter_option_parser_new(void *tool_arguments)
{
  GtSeqFilterArguments *arguments = tool_arguments;
  GtOption *option;
  GtOptionParser *op;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [sequence_file ...]",
                         "Filter the given sequence_file(s) and show the "
                         "results on stdout.");

  /* -minlength */
  option = gt_option_new_ulong("minlength",
                            "set minimum length a sequence must "
                            "have to pass the filter", &arguments->minlength,
                            GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, option);

  /* -maxlength */
  option = gt_option_new_ulong("maxlength", "set maximum length a sequence can "
                            "have to pass the filter", &arguments->maxlength,
                            GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_seqfilter_runner(int argc, const char **argv, int parsed_args,
                               void *tool_arguments, GtError *err)
{
  GtSeqFilterArguments *arguments = tool_arguments;
  GtBioseqIterator *bsi;
  GtBioseq *bioseq;
  unsigned long i;
  unsigned long long duplicates = 0, num_of_sequences = 0;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(tool_arguments);

  bsi = gt_bioseq_iterator_new(argc - parsed_args, argv + parsed_args);

  while (!(had_err = gt_bioseq_iterator_next(bsi, &bioseq, err)) && bioseq) {
    for (i = 0; i < gt_bioseq_number_of_sequences(bioseq); i++) {
      if ((arguments->minlength == GT_UNDEF_ULONG ||
           gt_bioseq_get_sequence_length(bioseq, i) >= arguments->minlength) &&
          (arguments->maxlength == GT_UNDEF_ULONG ||
           gt_bioseq_get_sequence_length(bioseq, i) <= arguments->maxlength)) {
        gt_fasta_show_entry(gt_bioseq_get_description(bioseq, i),
                            gt_bioseq_get_sequence(bioseq, i),
                            gt_bioseq_get_sequence_length(bioseq, i), 0);
      }
      else
        duplicates++;
      num_of_sequences++;
    }
    gt_bioseq_delete(bioseq);
  }

  /* show statistics */
  if (!had_err) {
    fprintf(stderr, "# %llu out of %llu sequences have been removed (%.3f%%)\n",
            duplicates, num_of_sequences,
            ((double) duplicates / num_of_sequences) * 100.0);
  }

  gt_bioseq_iterator_delete(bsi);

  return had_err;
}

GtTool* gt_seqfilter(void)
{
  return gt_tool_new(gt_seqfilter_arguments_new,
                  gt_seqfilter_arguments_delete,
                  gt_seqfilter_option_parser_new,
                  NULL,
                  gt_seqfilter_runner);
}
