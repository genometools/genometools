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
#include "libgtcore/alpha.h"
#include "libgtcore/error.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/score_matrix.h"
#include "libgtcore/unused.h"
#include "libgtext/blast_env.h"
#include "tools/gt_blastenv.h"

typedef struct {
  unsigned long q,
                k;
} ScorefastaArguments;

static void* gt_blastenv_arguments_new(void)
{
  return ma_malloc(sizeof (ScorefastaArguments));
}

static void gt_blastenv_arguments_delete(void *tool_arguments)
{
  ScorefastaArguments *arguments = tool_arguments;
  if (!arguments) return;
  ma_free(arguments);
}

static OptionParser* gt_blastenv_option_parser_new(void *tool_arguments)
{
  ScorefastaArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *o;
  assert(arguments);
  op = option_parser_new("[option ...] scorematrix_file w",
                         "Show the BlastP environment for sequence w (using "
                         "the given scormatrix_file).");
  o = option_new_ulong_min("q", "set q-gram length", &arguments->q, 4, 1);
  option_parser_add_option(op, o);
  o = option_new_ulong_min("k", "set minimum score", &arguments->k, 3, 1);
  option_parser_add_option(op, o);
  option_parser_set_min_max_args(op, 2, 2);
  return op;
}

static int gt_blastenv_runner(UNUSED int argc, const char **argv,
                              int parsed_args, void *tool_arguments,
                              Error *err)
{
  ScorefastaArguments *arguments = tool_arguments;
  ScoreMatrix *score_matrix;
  BlastEnv *blast_env = NULL;
  unsigned long wlen;
  char *w = NULL;
  Alpha *alpha = NULL;
  int had_err = 0;

  error_check(err);
  assert(arguments);

  if (!(score_matrix = score_matrix_new_read_protein(argv[parsed_args], err)))
    had_err = -1;

  if (!had_err) {
    /* store query sequence w */
    wlen = strlen(argv[parsed_args+1]);
    w = ma_malloc(wlen+1);
    strcpy(w, argv[parsed_args+1]);

    /* assign protein alphabet */
    alpha = alpha_new_protein();

    /* transform w according to protein alphabet */
    alpha_encode_seq(alpha, w, w, wlen);

    /* construct and show BlastP environment */
    blast_env = blast_env_new(w, wlen, alpha, arguments->q, arguments->k,
                              score_matrix);
    blast_env_show(blast_env);
  }

  /* free space */
  blast_env_delete(blast_env);
  score_matrix_delete(score_matrix);
  alpha_delete(alpha);
  ma_free(w);

  return had_err;
}

Tool* gt_blastenv(void)
{
  return tool_new(gt_blastenv_arguments_new,
                  gt_blastenv_arguments_delete,
                  gt_blastenv_option_parser_new,
                  NULL,
                  gt_blastenv_runner);
}
