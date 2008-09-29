/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "core/assert_api.h"
#include "core/ma.h"
#include "extended/qgramdist.h"
#include "extended/qgram.h"

unsigned long gt_calc_qgramdist(GtSeq *seq_a, GtSeq *seq_b, unsigned int q)
{
  unsigned long i, alphasize_to_the_power_of_q, *seq_a_profile, *seq_b_profile,
                dist = 0;
  const GtAlpha *gt_alpha_a, *gt_alpha_b;
  GtArray *seq_a_qgrams, *seq_b_qgrams;

  gt_assert(seq_a && seq_b);
  gt_alpha_a = gt_seq_get_alpha(seq_a);
  gt_alpha_b = gt_seq_get_alpha(seq_b);
  gt_assert(gt_alpha_is_compatible_with_alpha(gt_alpha_a, gt_alpha_b));
  alphasize_to_the_power_of_q = pow(gt_alpha_size(gt_alpha_a), q);

  seq_a_profile = gt_calloc(alphasize_to_the_power_of_q,
                            sizeof (unsigned long));
  seq_b_profile = gt_calloc(alphasize_to_the_power_of_q,
                            sizeof (unsigned long));

  seq_a_qgrams = gt_array_new(sizeof (unsigned long));
  seq_b_qgrams = gt_array_new(sizeof (unsigned long));

  gt_qgram_compute(seq_a_qgrams, gt_seq_get_encoded(seq_a),
                   gt_seq_length(seq_a), gt_alpha_size(gt_alpha_a), q);
  gt_assert(gt_array_size(seq_a_qgrams) == gt_seq_length(seq_a) - q + 1);
  gt_qgram_compute(seq_b_qgrams, gt_seq_get_encoded(seq_b),
                   gt_seq_length(seq_b), gt_alpha_size(gt_alpha_b), q);
  gt_assert(gt_array_size(seq_b_qgrams) == gt_seq_length(seq_b) - q + 1);

  for (i = 0; i < gt_array_size(seq_a_qgrams); i++)
    seq_a_profile[*(unsigned long*) gt_array_get(seq_a_qgrams, i)]++;
  for (i = 0; i < gt_array_size(seq_b_qgrams); i++)
    seq_b_profile[*(unsigned long*) gt_array_get(seq_b_qgrams, i)]++;

  /* compute distance */
  for (i = 0; i < alphasize_to_the_power_of_q; i++) {
    if (seq_a_profile[i] > seq_b_profile[i])
      dist += seq_a_profile[i] - seq_b_profile[i];
    else /* seq_a_profile[i] <= seq_b_profile[i] */
      dist += seq_b_profile[i] - seq_a_profile[i];
  }

  gt_array_delete(seq_b_qgrams);
  gt_array_delete(seq_a_qgrams);
  gt_free(seq_b_profile);
  gt_free(seq_a_profile);

  return dist;
}
