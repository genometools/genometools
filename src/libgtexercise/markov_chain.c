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

#include <assert.h>
#include <math.h>
#include <string.h>
#include "libgtcore/alpha.h"
#include "libgtcore/array2dim.h"
#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtexercise/markov_chain.h"

#define MINUSINFINITY   -99999.0

struct MarkovChain {
  Alpha *alpha;
  unsigned long num_of_states;
  double **transition_prob; /* log values */
};

MarkovChain* markov_chain_new(const char *states)
{
  unsigned long i, j;
  char characters[2];
  MarkovChain *mc;
  assert(states && strlen(states));
  /* alloc */
  mc = ma_malloc(sizeof *mc);
  mc->alpha = alpha_new();
  mc->num_of_states = strlen(states);
  array2dim_malloc(mc->transition_prob, mc->num_of_states, mc->num_of_states);
  /* set alphabet */
  characters[1] = '\0';
  for (i = 0 ; i < mc->num_of_states; i++) {
    characters[0] = states[i];
    alpha_add_mapping(mc->alpha, characters);
  }
  /* init */
  for (i = 0; i < mc->num_of_states; i++) {
    for (j = 0; j < mc->num_of_states; j++)
      markov_chain_set_transition_prob(mc, i, j, 1.0 / mc->num_of_states);
  }
  assert(markov_chain_is_valid(mc));
  return mc;
}

void markov_chain_delete(MarkovChain *mc)
{
  if (!mc) return;
  array2dim_delete(mc->transition_prob);
  alpha_delete(mc->alpha);
  ma_free(mc);
}

void markov_chain_set_transition_prob(MarkovChain *mc,
                                      unsigned long from_state_num,
                                      unsigned long to_state_num,
                                      double probability)
{
  assert(mc);
  assert(probability >= 0.0 && probability <= 1.0);
  assert(from_state_num < mc->num_of_states);
  assert(to_state_num < mc->num_of_states);
  if (probability == 0.0)
    mc->transition_prob[from_state_num][to_state_num] = MINUSINFINITY;
  else
    mc->transition_prob[from_state_num][to_state_num] = log(probability);
}

double markov_chain_get_transition_prob(const MarkovChain *mc,
                                        unsigned long from_state_num,
                                        unsigned long to_state_num)
{
  assert(mc);
  assert(from_state_num < mc->num_of_states);
  assert(to_state_num < mc->num_of_states);
  if (mc->transition_prob[from_state_num][to_state_num] == MINUSINFINITY)
    return 0.0;
  return  exp(mc->transition_prob[from_state_num][to_state_num]);
}

bool markov_chain_is_valid(const MarkovChain *mc)
{
  unsigned long i, j;
  assert(mc);
  for (i = 0; i < mc->num_of_states; i++) {
    double sum_of_probabilities = 0.0;
    for (j = 0; j < mc->num_of_states; j++)
      sum_of_probabilities += markov_chain_get_transition_prob(mc, i, j);
    if (!double_equals_one(sum_of_probabilities))
      return false;
  }
  return true;
}

int markov_chain_compute_prob(const MarkovChain *mc, double *prob,
                              const char *sequence, unsigned long seqlen,
                              Error *err)
{
  unsigned long i;
  unsigned int last_code = 0, cur_code;
  double logP = 0.0;
  int had_err = 0;
  error_check(err);
  assert(mc && sequence);
  assert(markov_chain_is_valid(mc));
  if (!alpha_char_is_valid(mc->alpha, sequence[0])) {
    error_set(err, "'%c' is not valid sequence character", sequence[0]);
    had_err = -1;
  }
  if (!had_err)
    last_code = alpha_encode(mc->alpha, sequence[0]);
  for (i = 1; !had_err && i < seqlen; i++) {
    if (!alpha_char_is_valid(mc->alpha, sequence[i])) {
      error_set(err, "'%c' is not valid sequence character", sequence[i]);
      had_err = -1;
      break;
    }
    cur_code = alpha_encode(mc->alpha, sequence[i]);
    if (mc->transition_prob[last_code][cur_code] == MINUSINFINITY) {
      logP = MINUSINFINITY;
      break;
    }
    else
      logP += mc->transition_prob[last_code][cur_code];
    last_code = cur_code;
  }
  if (!had_err) {
    if (logP == MINUSINFINITY)
      *prob = 0.0;
    else
      *prob = exp(logP);
  }
  return had_err;
}
