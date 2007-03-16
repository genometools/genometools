/*
  Copyright (c) 2005-2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <libgtext/coin_hmm.h>

HMM* coin_hmm_loaded(Env *env)
{
  HMM *hmm;

  /* create the HMM */
  hmm = hmm_new(COIN_NUM_OF_STATES, COIN_NUM_OF_SYMBOLS, env);

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

HMM* coin_hmm_fair(Env *env)
{
  HMM *hmm;

  /* create the HMM */
  hmm = hmm_new(COIN_NUM_OF_STATES, COIN_NUM_OF_SYMBOLS, env);

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

Alpha* coin_hmm_alpha(Env *env)
{
  Alpha *a = alpha_new(env);
  alpha_add_mapping(a, "Hh");
  alpha_add_mapping(a, "Tt");
  assert(alpha_size(a) == 2);
  return a;
}
