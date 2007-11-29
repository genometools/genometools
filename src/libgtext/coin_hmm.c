/*
  Copyright (c) 2005-2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2006 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/coin_hmm.h"

HMM* coin_hmm_loaded(void)
{
  HMM *hmm;

  /* create the HMM */
  hmm = hmm_new(COIN_NUM_OF_STATES, COIN_NUM_OF_SYMBOLS);

  /* set emmission probabilities */
  hmm_set_emission_probability(hmm, COIN_FAIR, HEAD, 0.5);
  hmm_set_emission_probability(hmm, COIN_FAIR, TAIL, 0.5);
  hmm_set_emission_probability(hmm, COIN_LOADED, HEAD, 0.75);
  hmm_set_emission_probability(hmm, COIN_LOADED, TAIL, 0.25);

  /* set transition probabilities */
  hmm_set_transition_probability(hmm, COIN_FAIR, COIN_LOADED, 0.1);
  hmm_set_transition_probability(hmm, COIN_LOADED, COIN_FAIR, 0.1);
  hmm_set_missing_transition_probabilities(hmm);
  assert(hmm_is_valid(hmm));

  return hmm;
}

HMM* coin_hmm_fair(void)
{
  HMM *hmm;

  /* create the HMM */
  hmm = hmm_new(COIN_NUM_OF_STATES, COIN_NUM_OF_SYMBOLS);

  /* set emmission probabilities */
  hmm_set_emission_probability(hmm, COIN_FAIR, HEAD, 0.5);
  hmm_set_emission_probability(hmm, COIN_FAIR, TAIL, 0.5);
  hmm_set_emission_probability(hmm, COIN_LOADED, HEAD, 0.5);
  hmm_set_emission_probability(hmm, COIN_LOADED, TAIL, 0.5);

  /* set transition probabilities */
  hmm_set_transition_probability(hmm, COIN_FAIR, COIN_LOADED, 0.5);
  hmm_set_transition_probability(hmm, COIN_LOADED, COIN_FAIR, 0.5);
  hmm_set_missing_transition_probabilities(hmm);
  assert(hmm_is_valid(hmm));

  return hmm;
}

Alpha* coin_hmm_alpha(void)
{
  Alpha *a = alpha_new();
  alpha_add_mapping(a, "Hh");
  alpha_add_mapping(a, "Tt");
  assert(alpha_size(a) == 2);
  return a;
}
