/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <math.h>
#include "qgram.h"

unsigned long qgram_encode(const char *w, unsigned long q,
                           unsigned long alphabet_size)
{
  unsigned long i, qgram_code;
  assert(w);
  qgram_code = w[0];
  for (i = 1; i < q; i++) {
    if (alphabet_size == 4)
      qgram_code = qgram_code << 2 | w[i];
    else
      qgram_code = qgram_code * alphabet_size + w[i];
  }
  return qgram_code;
}

unsigned long qgram_step(unsigned long current_code, char previous, char next,
                         unsigned long alphabet_size,
                         unsigned long
                         alpha_size_raised_to_the_power_of_q_minus_1)

{
  unsigned long next_code;
  next_code = (current_code - previous *
               alpha_size_raised_to_the_power_of_q_minus_1) * alphabet_size +
              next;
  return next_code;
}

void qgram_compute(Array *qgrams, const char *encoded_seq,
                   unsigned long seqlen, unsigned long alpha_size,
                   unsigned int q, Env *env)
{
  unsigned long i, code, alpha_size_raised_to_the_power_of_q_minus_1;
  assert(qgrams && encoded_seq && alpha_size && q);
  if (seqlen >= q) {
    alpha_size_raised_to_the_power_of_q_minus_1 = pow(alpha_size, q-1);
    code = qgram_encode(encoded_seq, q, alpha_size);
    array_add(qgrams, code, env);
    i = 0;
    while (i + q < seqlen) {
      code = qgram_step(code, encoded_seq[i], encoded_seq[i+q], alpha_size,
                        alpha_size_raised_to_the_power_of_q_minus_1);
      array_add(qgrams, code, env);
      i++;
    }
  }
}
