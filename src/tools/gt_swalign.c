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

#include "core/bioseq.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/score_function.h"
#include "core/unused_api.h"
#include "core/xansi.h"
#include "extended/alignment.h"
#include "extended/swalign.h"
#include "tools/gt_swalign.h"

#define DEFAULT_INDELSCORE  -3

typedef struct {
  int indelscore;
}  SWAlignArguments;

static void* gt_swalign_arguments_new(void)
{
  return gt_calloc(1, sizeof (SWAlignArguments));
}

static void gt_swalign_arguments_delete(void *tool_arguments)
{
  SWAlignArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static OptionParser* gt_swalign_opion_parser_new(void *tool_arguments)
{
  SWAlignArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *o;
  assert(arguments);
  op = option_parser_new("[option ...] scorematrix seq_file_1 seq_file_2",
                         "Locally align each sequence in seq_file_1 "
                         "with each sequence in seq_file_2.");
  o = option_new_int("indelscore", "set the score used for "
                     "insertions/deletions", &arguments->indelscore,
                     DEFAULT_INDELSCORE);
  option_parser_add_option(op, o);
  option_parser_set_min_max_args(op, 3, 3);
  return op;
}

static int gt_swalign_runner(GT_UNUSED int argc, const char **argv,
                             int parsed_args, void *tool_arguments,
                             GtError *err)
{
  SWAlignArguments *arguments = tool_arguments;
  GT_Bioseq *gt_bioseq_1 = NULL, *gt_bioseq_2 = NULL;
  GT_ScoreFunction *score_function = NULL;
  GT_ScoreMatrix *scorematrix;
  unsigned long i, j;
  int had_err = 0;
  Alignment *a;
  gt_error_check(err);
  assert(arguments);

  /* init */
  /* XXX: make this more flexible */
  scorematrix  = gt_score_matrix_new_read_protein(argv[parsed_args], err);
  if (scorematrix) {
    score_function = gt_score_function_new(scorematrix, arguments->indelscore,
                                                   arguments->indelscore);
    gt_bioseq_1 = gt_bioseq_new(argv[parsed_args+1], err);
    if (!gt_bioseq_1)
      had_err = -1;
    if (!had_err) {
      gt_bioseq_2 = gt_bioseq_new(argv[parsed_args+2], err);
      if (!gt_bioseq_2)
        had_err = -1;
    }

    if (!had_err) {
      /* aligning all sequence combinations */
      for (i = 0; i < gt_bioseq_number_of_sequences(gt_bioseq_1); i++) {
        for (j = 0; j < gt_bioseq_number_of_sequences(gt_bioseq_2); j++) {
          a = swalign(gt_bioseq_get_seq(gt_bioseq_1, i),
                      gt_bioseq_get_seq(gt_bioseq_2, j), score_function);
          if (a) {
            alignment_show(a, stdout);
            xputchar('\n');
            alignment_delete(a);
          }
        }
      }
    }
  }

  /* free */
  gt_bioseq_delete(gt_bioseq_2);
  gt_bioseq_delete(gt_bioseq_1);
  gt_score_function_delete(score_function);

  return had_err;
}

Tool* gt_swalign(void)
{
  return tool_new(gt_swalign_arguments_new,
                  gt_swalign_arguments_delete,
                  gt_swalign_opion_parser_new,
                  NULL,
                  gt_swalign_runner);
}
