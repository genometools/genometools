/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "extended/qgram.h"

unsigned long gt_qgram_encode(const char *w, unsigned long q,
                           unsigned long alphabet_size)
{
  unsigned long i, qgram_code;
  gt_assert(w);
  qgram_code = w[0];
  for (i = 1; i < q; i++) {
    if (alphabet_size == 4)
      qgram_code = qgram_code << 2 | w[i];
    else
      qgram_code = qgram_code * alphabet_size + w[i];
  }
  return qgram_code;
}

unsigned long gt_qgram_step(unsigned long current_code, char previous,
                            char next, unsigned long alphabet_size,
                            unsigned long
                            alpha_size_raised_to_the_power_of_q_minus_1)

{
  unsigned long next_code;
  next_code = (current_code - previous *
               alpha_size_raised_to_the_power_of_q_minus_1) * alphabet_size +
               next;
  return next_code;
}

void gt_qgram_compute(GtArray *qgrams, const char *encoded_seq,
                   unsigned long seqlen, unsigned long gt_alpha_size,
                   unsigned int q)
{
  unsigned long i, code, alpha_size_raised_to_the_power_of_q_minus_1;
  gt_assert(qgrams && encoded_seq && gt_alpha_size && q);
  if (seqlen >= q) {
    alpha_size_raised_to_the_power_of_q_minus_1 = pow(gt_alpha_size, q-1);
    code = gt_qgram_encode(encoded_seq, q, gt_alpha_size);
    gt_array_add(qgrams, code);
    i = 0;
    while (i + q < seqlen) {
      code = gt_qgram_step(code, encoded_seq[i], encoded_seq[i+q],
                           gt_alpha_size,
                           alpha_size_raised_to_the_power_of_q_minus_1);
      gt_array_add(qgrams, code);
      i++;
    }
  }
}

void gt_qgram_decode(char *qgram, unsigned long code, unsigned long q,
                  GtAlpha *alpha)
{
  unsigned int alphabet_size, c = 0;
  unsigned long i;
  gt_assert(qgram && q && alpha);
  alphabet_size = gt_alpha_size(alpha);
  for (i = q; i > 0; i--) {
    c = code % alphabet_size;
    code = (code - c) / alphabet_size;
    qgram[i-1] = gt_alpha_decode(alpha, c);
  }
}
