/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef HMM_H
#define HMM_H

#include "core/error.h"

/* The Hidden Markov Model (HMM) class */
typedef struct GtHMM GtHMM;

GtHMM*   gt_hmm_new(unsigned int num_of_states, unsigned int num_of_symbols);
void   gt_hmm_set_initial_state_probability(GtHMM*, unsigned int state_num,
                                         double probability);
double gt_hmm_get_initial_state_probability(const GtHMM*,
                                            unsigned int state_num);
void   gt_hmm_set_transition_probability(GtHMM*, unsigned int from_state_num,
                                         unsigned int to_state_num,
                                         double probability);
void   gt_hmm_set_missing_transition_probabilities(GtHMM*);
double gt_hmm_get_transition_probability(const GtHMM*,
                                         unsigned int from_state_num,
                                         unsigned int to_state_num);
void   gt_hmm_set_emission_probability(GtHMM*, unsigned int state_num,
                                       unsigned int symbol_num,
                                       double probability);
double gt_hmm_get_emission_probability(const GtHMM*, unsigned int state_num,
                                                unsigned int symbol_num);
/* initialize the GtHMM with completly random values */
void   gt_hmm_init_random(GtHMM*);
/* Viterbi algorithm */
void   gt_hmm_decode(const GtHMM*, unsigned int *state_sequence,
                  const unsigned int *emissions, unsigned int num_of_emissions);
/* Forward algorithm, returns log(P(emissions)) */
double gt_hmm_forward(const GtHMM*, const unsigned int *emissions,
                   unsigned int num_of_emissions);
/* Backward algorithm, returns log(P(emissions)) */
double gt_hmm_backward(const GtHMM*, const unsigned int *emissions,
                    unsigned int num_of_emissions);
void   gt_hmm_emit(GtHMM*, unsigned long num_of_emissions,
                void (*proc_emission)(unsigned int symbol, void *data),
                void *data);
bool   gt_hmm_is_valid(const GtHMM*);
/* returns the RMSD of two GtHMMs */
double gt_hmm_rmsd(const GtHMM*, const GtHMM*);
void   gt_hmm_show(const GtHMM*, FILE*);
int    gt_hmm_unit_test(GtError*);
void   gt_hmm_delete(GtHMM*);

#endif
