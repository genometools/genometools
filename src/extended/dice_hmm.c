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

#include "core/assert_api.h"
#include "extended/dice_hmm.h"

GtHMM* gt_dice_hmm_loaded(void)
{
  GtHMM *hmm;

  /* create the GtHMM */
  hmm = gt_hmm_new(DICE_NUM_OF_STATES, DICE_NUM_OF_SYMBOLS);

  /* set emmission probabilities */
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, ONE,   1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, TWO,   1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, THREE, 1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, FOUR,  1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, FIVE,  1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, SIX,   1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, ONE,   1.0/10);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, TWO,   1.0/10);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, THREE, 1.0/10);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, FOUR,  1.0/10);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, FIVE,  1.0/10);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, SIX,   1.0/2);

  /* set transition probabilities */
  gt_hmm_set_transition_probability(hmm, DICE_FAIR, DICE_LOADED, 0.05);
  gt_hmm_set_transition_probability(hmm, DICE_LOADED, DICE_FAIR, 0.1);
  gt_hmm_set_transition_probability(hmm, DICE_FAIR, DICE_FAIR, 0.95);
  gt_hmm_set_transition_probability(hmm, DICE_LOADED, DICE_LOADED, 0.9);
  gt_assert(gt_hmm_is_valid(hmm));

  return hmm;
}

GtHMM* gt_dice_hmm_fair(void)
{
  GtHMM *hmm;

  /* create the GtHMM */
  hmm = gt_hmm_new(DICE_NUM_OF_STATES, DICE_NUM_OF_SYMBOLS);

  /* set emmission probabilities */
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, ONE,   1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, TWO,   1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, THREE, 1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, FOUR,  1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, FIVE,  1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_FAIR, SIX,   1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, ONE,   1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, TWO,   1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, THREE, 1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, FOUR,  1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, FIVE,  1.0/6);
  gt_hmm_set_emission_probability(hmm, DICE_LOADED, SIX,   1.0/6);

  /* set transition probabilities */
  gt_hmm_set_transition_probability(hmm, DICE_FAIR, DICE_LOADED, 0.5);
  gt_hmm_set_transition_probability(hmm, DICE_LOADED, DICE_FAIR, 0.5);
  gt_hmm_set_transition_probability(hmm, DICE_FAIR, DICE_FAIR, 0.5);
  gt_hmm_set_transition_probability(hmm, DICE_LOADED, DICE_LOADED, 0.5);
  gt_assert(gt_hmm_is_valid(hmm));

  return hmm;
}

GtAlphabet* gt_dice_hmm_alphabet(void)
{
  GtAlphabet *a = gt_alphabet_new_empty();
  gt_alphabet_add_mapping(a, "1");
  gt_alphabet_add_mapping(a, "2");
  gt_alphabet_add_mapping(a, "3");
  gt_alphabet_add_mapping(a, "4");
  gt_alphabet_add_mapping(a, "5");
  gt_alphabet_add_mapping(a, "6");
  gt_assert(gt_alphabet_size(a) == 6);
  return a;
}
