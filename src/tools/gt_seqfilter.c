/*
  Copyright (c) 2008-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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
#include "core/output_file_api.h"
#include "core/option_api.h"
#include "core/undef_api.h"
#include "tools/gt_seqfilter.h"

typedef struct {
  unsigned long minlength,
                maxlength,
                maxseqnum,
                width;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} SeqFilterArguments;

static void* gt_seqfilter_arguments_new(void)
{
  SeqFilterArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_seqfilter_arguments_delete(void *tool_arguments)
{
  SeqFilterArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_seqfilter_option_parser_new(void *tool_arguments)
{
  SeqFilterArguments *arguments = tool_arguments;
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

  /* -maxseqnum */
  option = gt_option_new_ulong("maxseqnum", "set the maximum number of "
                               "sequences which can pass the filter",
                               &arguments->maxseqnum, GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, option);

  /* -width */
  option = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, option);

  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  return op;
}

static int gt_seqfilter_runner(int argc, const char **argv, int parsed_args,
                               void *tool_arguments, GtError *err)
{
  SeqFilterArguments *arguments = tool_arguments;
  GtBioseqIterator *bsi;
  GtBioseq *bioseq;
  unsigned long i;
  unsigned long long passed = 0, filtered = 0, num_of_sequences = 0;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(tool_arguments);

  bsi = gt_bioseq_iterator_new(argc - parsed_args, argv + parsed_args);

  while (!(had_err = gt_bioseq_iterator_next(bsi, &bioseq, err)) && bioseq) {
    for (i = 0; i < gt_bioseq_number_of_sequences(bioseq); i++) {
      char *seq;
      if ((arguments->minlength == GT_UNDEF_ULONG ||
           gt_bioseq_get_sequence_length(bioseq, i) >= arguments->minlength) &&
          (arguments->maxlength == GT_UNDEF_ULONG ||
           gt_bioseq_get_sequence_length(bioseq, i) <= arguments->maxlength) &&
          (arguments->maxseqnum == GT_UNDEF_ULONG ||
           passed + 1 <= arguments->maxseqnum)) {
        seq = gt_bioseq_get_sequence(bioseq, i);
        gt_fasta_show_entry(gt_bioseq_get_description(bioseq, i),
                            seq,
                            gt_bioseq_get_sequence_length(bioseq, i),
                            arguments->width, arguments->outfp);
        gt_free(seq);
        passed++;
      }
      else
        filtered++;
      num_of_sequences++;
    }
    gt_bioseq_delete(bioseq);
  }

  /* show statistics */
  if (!had_err) {
    gt_assert(passed + filtered == num_of_sequences);
    fprintf(stderr, "# %llu out of %llu sequences have been removed (%.3f%%)\n",
            filtered, num_of_sequences,
            ((double) filtered / num_of_sequences) * 100.0);
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
