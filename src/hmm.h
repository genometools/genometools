/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef HMM_H
#define HMM_H

/* The Hidden Markov Model (HMM) class */
typedef struct HMM HMM;

HMM*         hmm_new(unsigned int num_of_states, unsigned int num_of_symbols);
void         hmm_set_initial_state_probability(HMM*, unsigned int state_num,
                                               double probability);
double       hmm_get_initial_state_probability(const HMM*,
                                               unsigned int state_num);
void         hmm_set_transition_probability(HMM*,
                                            unsigned int from_state_num,
                                            unsigned int to_state_num,
                                            double probability);
void         hmm_set_missing_transition_probabilities(HMM*);
double       hmm_get_transition_probability(const HMM*,
                                            unsigned int from_state_num,
                                            unsigned int to_state_num);
void         hmm_set_emission_probability(HMM*, unsigned int state_num,
                                          unsigned int symbol_num,
                                          double probability);
double       hmm_get_emission_probability(const HMM*, unsigned int state_num,
                                          unsigned int symbol_num);
/* initialize the HMM with completly random values */
void         hmm_init_random(HMM*);
/* Viterbi algorithm */
void         hmm_decode(const HMM*, unsigned int *state_sequence,
                        const unsigned int *emissions,
                        unsigned int num_of_emissions);
/* Forward algorithm, returns log(P(emissions)) */
double       hmm_forward(const HMM*, const unsigned int *emissions,
                         unsigned int num_of_emissions);
/* Backward algorithm, returns log(P(emissions)) */
double       hmm_backward(const HMM*, const unsigned int *emissions,
                          unsigned int num_of_emissions);
void         hmm_emit(HMM*, unsigned long num_of_emissions,
                      void (*proc_emission)(unsigned int symbol, void *data),
                      void *data);
unsigned int hmm_is_valid(const HMM*);
/* returns the RMSD of two HMMs */
double       hmm_rmsd(const HMM*, const HMM*);
void         hmm_show(const HMM*, FILE*);
int          hmm_unit_test(void);
void         hmm_free(HMM*);

#endif
