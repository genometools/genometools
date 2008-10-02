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

#include <math.h>
#include <string.h>
#include "core/unused_api.h"
#include "exercise/markov_chain_parsing.h"
#include "tools/gt_markovchain.h"

static GtOptionParser* gt_markovchain_option_parser_new(GT_UNUSED
                                                      void *tool_arguments)
{
  GtOptionParser *op;
  op = gt_option_parser_new("[option ...] markov_chain_file sequence",
                         "Compute the probability of sequence given the markov "
                         "chain in markov_chain_file.");
  gt_option_parser_set_min_max_args(op, 2, 2);
  return op;
}

static int gt_markovchain_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args,
                                 GT_UNUSED void *tool_arguments, GtError *err)
{
  GtMarkovChain *mc;
  double P;
  int had_err = 0;

  gt_error_check(err);

  if (!(mc = gt_markov_chain_parse(argv[parsed_args], err)))
    had_err = -1;

  if (!had_err) {
    unsigned long seqlen = strlen(argv[parsed_args+1]);
    gt_assert(gt_markov_chain_is_valid(mc));
    had_err = gt_markov_chain_compute_prob(mc, &P, argv[parsed_args+1], seqlen,
                                        err);
  }

  if (!had_err)
    printf("P=%.10f\n", P);

  gt_markov_chain_delete(mc);

  return had_err;
}

GtTool* gt_markovchain(void)
{
  return gt_tool_new(NULL,
                  NULL,
                  gt_markovchain_option_parser_new,
                  NULL,
                  gt_markovchain_runner);
}
