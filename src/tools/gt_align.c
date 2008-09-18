/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/bioseq.h"
#include "core/option.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/xansi.h"
#include "extended/align.h"
#include "extended/alignment.h"
#include "tools/gt_align.h"

typedef struct {
  bool all;
} AlignArguments;

static void* gt_align_arguments_new(void)
{
  return gt_calloc(1, sizeof (AlignArguments));
}

static void gt_align_arguments_delete(void *tool_arguments)
{
  AlignArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static GtOptionParser* gt_align_option_parser_new(void *tool_arguments)
{
  AlignArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  assert(arguments);

  op = gt_option_parser_new("[option ...] seq_file_1 seq_file_2",
                            "Globally align each sequence in seq_file_1 with "
                            "each sequence in seq_file_2.");
  option = gt_option_new_bool("all",
                              "show all optimal alignments instead of just "
                              "one", &arguments->all, false);
  gt_option_parser_add_option(op, option);
  gt_option_parser_set_min_max_args(op, 2, 2);
  return op;
}

static void show_alignment(const GtAlignment *a, GT_UNUSED void *data)
{
  assert(a && !data);
  gt_alignment_show(a, stdout);
  xputchar('\n');
}

static void show_aligns(unsigned long aligns, GT_UNUSED void *data)
{
  assert(aligns && !data);
  printf("number of optimal alignments: %lu\n\n", aligns);
}

static int gt_align_runner(GT_UNUSED int argc, const char **argv,
                           int parsed_args, void *tool_arguments, GtError *err)
{
  AlignArguments *arguments = tool_arguments;
  GtBioseq *gt_bioseq_1, *gt_bioseq_2 = NULL;
  unsigned long i, j;
  int had_err = 0;
  GtAlignment *a;
  gt_error_check(err);
  assert(arguments);

  /* init */
  gt_bioseq_1 = gt_bioseq_new(argv[parsed_args], err);
  if (!gt_bioseq_1)
    had_err = -1;
  if (!had_err) {
    gt_bioseq_2 = gt_bioseq_new(argv[parsed_args+1], err);
    if (!gt_bioseq_2)
      had_err = -1;
  }

  /* aligning all sequence combinations */
  if (!had_err) {
    for (i = 0; i < gt_bioseq_number_of_sequences(gt_bioseq_1); i++) {
      for (j = 0; j < gt_bioseq_number_of_sequences(gt_bioseq_2); j++) {
        if (arguments->all) {
          align_all(gt_bioseq_get_sequence(gt_bioseq_1, i),
                    gt_bioseq_get_sequence_length(gt_bioseq_1, i),
                    gt_bioseq_get_sequence(gt_bioseq_2, j),
                    gt_bioseq_get_sequence_length(gt_bioseq_2, j),
                    show_alignment, show_aligns, NULL);
        }
        else {
          a = align(gt_bioseq_get_sequence(gt_bioseq_1, i),
                    gt_bioseq_get_sequence_length(gt_bioseq_1, i),
                    gt_bioseq_get_sequence(gt_bioseq_2, j),
                    gt_bioseq_get_sequence_length(gt_bioseq_2, j));
          gt_alignment_show(a, stdout);
          xputchar('\n');
          gt_alignment_delete(a);
        }
      }
    }
  }

  /* free */
  gt_bioseq_delete(gt_bioseq_2);
  gt_bioseq_delete(gt_bioseq_1);

  return had_err;
}

GtTool* gt_align(void)
{
  return gt_tool_new(gt_align_arguments_new,
                  gt_align_arguments_delete,
                  gt_align_option_parser_new,
                  NULL,
                  gt_align_runner);
}
