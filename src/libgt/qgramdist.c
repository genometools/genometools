/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <math.h>
#include "array.h"
#include "qgramdist.h"
#include "qgram.h"
#include "xansi.h"

unsigned long qgramdist(Seq *seq_a, Seq *seq_b, unsigned int q)
{
  unsigned long i, alphasize_to_the_power_of_q, *seq_a_profile, *seq_b_profile,
                dist = 0;
  const Alpha *alpha_a, *alpha_b;
  Array *seq_a_qgrams, *seq_b_qgrams;

  assert(seq_a && seq_b);
  alpha_a = seq_get_alpha(seq_a);
  alpha_b = seq_get_alpha(seq_b);
  assert(alpha_is_compatible_with_alpha(alpha_a, alpha_b));
  alphasize_to_the_power_of_q = pow(alpha_size(alpha_a), q);

  seq_a_profile = xcalloc(alphasize_to_the_power_of_q, sizeof (unsigned long));
  seq_b_profile = xcalloc(alphasize_to_the_power_of_q, sizeof (unsigned long));

  seq_a_qgrams = array_new(sizeof (unsigned long));
  seq_b_qgrams = array_new(sizeof (unsigned long));

  qgram_compute(seq_a_qgrams, seq_get_encoded(seq_a), seq_length(seq_a),
                alpha_size(alpha_a), q);
  assert(array_size(seq_a_qgrams) == seq_length(seq_a) - q + 1);
  qgram_compute(seq_b_qgrams, seq_get_encoded(seq_b), seq_length(seq_b),
                alpha_size(alpha_b), q);
  assert(array_size(seq_b_qgrams) == seq_length(seq_b) - q + 1);

  for (i = 0; i < array_size(seq_a_qgrams); i++)
    seq_a_profile[*(unsigned long*) array_get(seq_a_qgrams, i)]++;
  for (i = 0; i < array_size(seq_b_qgrams); i++)
    seq_b_profile[*(unsigned long*) array_get(seq_b_qgrams, i)]++;

  /* compute distance */
  for (i = 0; i < alphasize_to_the_power_of_q; i++) {
    if (seq_a_profile[i] > seq_b_profile[i])
      dist += seq_a_profile[i] - seq_b_profile[i];
    else /* seq_a_profile[i] <= seq_b_profile[i] */
      dist += seq_b_profile[i] - seq_a_profile[i];
  }

  array_delete(seq_b_qgrams);
  array_delete(seq_a_qgrams);
  free(seq_b_profile);
  free(seq_a_profile);

  return dist;
}
