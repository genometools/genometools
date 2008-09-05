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

#ifndef MARKOV_CHAIN_H
#define MARKOV_CHAIN_H

#include <stdbool.h>

typedef struct MarkovChain MarkovChain;

MarkovChain* markov_chain_new(const char *states);
void         markov_chain_delete(MarkovChain*);
void         markov_chain_set_transition_prob(MarkovChain*,
                                              unsigned long from_state_num,
                                              unsigned long to_state_num,
                                              double probability);
double       markov_chain_get_transition_prob(const MarkovChain*,
                                              unsigned long from_state_num,
                                              unsigned long to_state_num);
bool         markov_chain_is_valid(const MarkovChain*);
int          markov_chain_compute_prob(const MarkovChain*, double *prob,
                                       const char *sequence,
                                       unsigned long seqlen, Error*);

#endif
