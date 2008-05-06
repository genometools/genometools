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

#include "libgtcore/bioseq.h"
#include "libgtcore/option.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"
#include "libgtext/align.h"
#include "libgtext/alignment.h"
#include "tools/gt_align.h"

typedef struct {
  bool all;
} AlignArguments;

static void* gt_align_arguments_new(void)
{
  return ma_calloc(1, sizeof (AlignArguments));
}

static void gt_align_arguments_delete(void *tool_arguments)
{
  AlignArguments *arguments = tool_arguments;
  if (!arguments) return;
  ma_free(arguments);
}

static OptionParser* gt_align_option_parser_new(void *tool_arguments)
{
  AlignArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *option;
  assert(arguments);

  op = option_parser_new("[option ...] seq_file_1 seq_file_2",
                         "Globally align each sequence in seq_file_1 with each "
                         "sequence in seq_file_2.");
  option = option_new_bool("all", "show all optimal alignments instead of just "
                           "one", &arguments->all, false);
  option_parser_add_option(op, option);
  option_parser_set_min_max_args(op, 2, 2);
  return op;
}

static void show_alignment(const Alignment *a, UNUSED void *data)
{
  assert(a && !data);
  alignment_show(a, stdout);
  xputchar('\n');
}

static void show_aligns(unsigned long aligns, UNUSED void *data)
{
  assert(aligns && !data);
  printf("number of optimal alignments: %lu\n\n", aligns);
}

static int gt_align_runner(UNUSED int argc, const char **argv, int parsed_args,
                           void *tool_arguments, Error *err)
{
  AlignArguments *arguments = tool_arguments;
  Bioseq *bioseq_1, *bioseq_2 = NULL;
  unsigned long i, j;
  int had_err = 0;
  Alignment *a;
  error_check(err);
  assert(arguments);

  /* init */
  bioseq_1 = bioseq_new(argv[parsed_args], err);
  if (!bioseq_1)
    had_err = -1;
  if (!had_err) {
    bioseq_2 = bioseq_new(argv[parsed_args+1], err);
    if (!bioseq_2)
      had_err = -1;
  }

  /* aligning all sequence combinations */
  if (!had_err) {
    for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
      for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
        if (arguments->all) {
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
          alignment_delete(a);
        }
      }
    }
  }

  /* free */
  bioseq_delete(bioseq_2);
  bioseq_delete(bioseq_1);

  return had_err;
}

Tool* gt_align(void)
{
  return tool_new(gt_align_arguments_new,
                  gt_align_arguments_delete,
                  gt_align_option_parser_new,
                  NULL,
                  gt_align_runner);
}
