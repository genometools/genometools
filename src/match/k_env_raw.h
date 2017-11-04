/*
  Copyright (c) 2017 Julian Elvers <julian.elvers@gmail.com>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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

#ifndef K_ENV_RAW_H
#define K_ENV_RAW_H

#include <stdbool.h>
#include "core/error_api.h"
#include "core/types_api.h"
#include "core/score_matrix.h"
#include "k_env_alphabet.h"

/* Struct describing an element of the k-environment of a specific kmer.
   <kmer_neighbor> is the encoded Kmer, which has a score of <score> with
   the reference kmer. */
typedef struct
{
  GtCodetype kmer_neighbor;
  int score;
} GtKenvRawElement;

/* Struct introducing the GtKenvGenerator class, used to build k_environments
   for given encoded kmers or sequences. */
typedef struct GtKenvGenerator GtKenvGenerator;

/* Returns a new empty <GtKenvGenerator> object giving the frame to
   build a k environment. */
GtKenvGenerator *gt_kenv_generator_new(const GtKenvAlphabet *kenv_alphabet,
                                       unsigned int q_val, int th_val,
                                       bool allow_x_val, bool preprocess,
                                       GtError *err);

/* Deletes the given <GtKenvGenerator>. */
void gt_kenv_generator_delete(GtKenvGenerator *kenv_gen);

/* Returns a pointer to the memory of <GtKenvRawElement> entries, representing
   the k_environment for the given <refcode>. */
const GtKenvRawElement *gt_kenv_raw_element_list(GtUword *length,
                                GtCodetype refcode, GtKenvGenerator *kenv_gen);

/* Returns the code for the given <word>, encoded by the alphabet for the
   <kenv_gen>*/
GtCodetype gt_kenv_get_kmer_code(const GtKenvGenerator *kenv_gen,
                                 const GtUchar *word);

/* Prepares the <kenv_gen> buffers for the given <kmer> to enable
   the usage of the further methods. */
void gt_kenv_generator_fill_buffers_from_sequence(GtKenvGenerator *kenv_gen,
                                                  const GtUchar *kmer);

/* Returns the value for <q> stored inside <kenv_gen>. <q> is the defined
   length of the kmeres. */
unsigned int gt_kenv_generator_get_q(const GtKenvGenerator *kenv_gen);

/* Returns the value for <th> stored inside <kenv_gen>. <th> is the defined
   score threshold of the kmeres inside the environment. */
int gt_kenv_generator_get_th(const GtKenvGenerator *kenv_gen);

/* Returns the number of bits required for storing a score value in the
   range from score_threshold to max_score_qgram */
size_t gt_kenv_generator_get_score_bits(const GtKenvGenerator *kenv_gen);

/* Returns the value for <allow_w> stored inside <kenv_gen>. This is the chosen
   handling for wildcard amino acid X. */
bool gt_kenv_generator_get_allow_x(const GtKenvGenerator *kenv_gen);

/* Returns the array containing the maximum score for each character of the
   Alphabet within <kenv_gen>. */
const int *gt_kenv_generator_get_max_score_each_char(
                                              const GtKenvGenerator *kenv_gen);

/* Returns the highest score possible for a given q-gram. */
int gt_kenv_generator_get_max_score_qgram(const GtKenvGenerator *kenv_gen);

/* return score for two encoded sequences of the given length with respect
   to the kenvironment */

int gt_kenv_eval_score(const GtKenvGenerator *kenv_gen,
                       const GtUchar *a_encoded,
                       const GtUchar *b_encoded,
                       GtUword length);

/* return reference to a score_matrix (which can be NULL) */

const GtScoreMatrix *gt_kenv_score_matrix(const GtKenvGenerator *kenv_gen);

#endif
