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
#include "core/alpha.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "extended/scorefasta.h"
#include "tools/gt_scorefasta.h"

typedef struct {
  unsigned long q;
} ScorefastaArguments;

static void* gt_scorefasta_arguments_new(void)
{
  return gt_malloc(sizeof (ScorefastaArguments));
}

static void gt_scorefasta_arguments_delete(void *tool_arguments)
{
  ScorefastaArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static OptionParser* gt_scorefasta_option_parser_new(void *tool_arguments)
{
  ScorefastaArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *o;
  assert(arguments);
  op = option_parser_new("[option ...] u w",
                         "Compute scorefasta for DNA sequences u and w).");
  o = option_new_ulong_min("q", "set q-gram length", &arguments->q, 3, 1);
  option_parser_add_option(op, o);
  option_parser_set_min_max_args(op, 2, 2);
  return op;
}

static int gt_scorefasta_runner(GT_UNUSED int argc, const char **argv,
                                int parsed_args, void *tool_arguments,
                                GT_UNUSED GtError *err)
{
  ScorefastaArguments *arguments = tool_arguments;
  unsigned long ulen, wlen;
  char *u, *w;
  GT_Alpha *alpha;

  gt_error_check(err);
  assert(arguments);

  /* store database sequence u and query sequence w */
  ulen = strlen(argv[parsed_args]);
  wlen = strlen(argv[parsed_args+1]);
  u = gt_malloc(ulen+1);
  w = gt_malloc(wlen+1);
  strcpy(u, argv[parsed_args]);
  strcpy(w, argv[parsed_args+1]);

  /* assign DNA alphabet */
  alpha = gt_alpha_new_dna();

  /* transform u and w according to DNA alphabet */
  gt_alpha_encode_seq(alpha, u, u, ulen);
  gt_alpha_encode_seq(alpha, w, w, wlen);

  printf("scorefasta=%lu\n", scorefasta(u, ulen, w, wlen, arguments->q,
                                        gt_alpha_size(alpha)));

  /* free space */
  gt_alpha_delete(alpha);
  gt_free(u);
  gt_free(w);

  return 0;
}

Tool* gt_scorefasta(void)
{
  return tool_new(gt_scorefasta_arguments_new,
                  gt_scorefasta_arguments_delete,
                  gt_scorefasta_option_parser_new,
                  NULL,
                  gt_scorefasta_runner);
}
