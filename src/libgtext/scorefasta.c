/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/array.h"
#include "libgtcore/ma.h"
#include "libgtext/qgram.h"

unsigned long scorefasta(const char *u, unsigned long ulen,
                         const char *w, unsigned long wlen,
                         unsigned long q, unsigned long alphabet_size)
{
  Array **h;
  unsigned long i, j, code, hsize, r_raised_to_the_power_of_q_minus_1,
                *count, maxcount = 0;

  /* some checks */
  if (q > ulen || q > ulen)
    return 0;

  /* preprocess function h */
  hsize = pow(alphabet_size, q);
  h = ma_malloc(sizeof(Array*) * hsize);
  for (i = 0; i < hsize; i++)
    h[i] = array_new(sizeof (unsigned long));

  /* compute code of first q-gram of w */
  code = qgram_encode(w, q, alphabet_size);

  /* store first encoding in h */
  i = 0;
  array_add(h[code], i);

  r_raised_to_the_power_of_q_minus_1 = pow(alphabet_size, q-1);
  for (i = 1; i < wlen - q + 1; i++) {
    /* update code */
    code = qgram_step(code, w[i-1], w[i+q-1], alphabet_size,
                      r_raised_to_the_power_of_q_minus_1);
    /* store startposition in h[code] */
    array_add(h[code], i);
  }

  /* final phase */

  /* init count */
  count = ma_calloc((ulen + wlen + 1), sizeof *count);

  /* compute code of first q-gram in u */
  code = qgram_encode(u, q, alphabet_size);

  /* increase count */
  for (i = 0; i < array_size(h[code]); i++)
    count[wlen - 1 - *(unsigned long*) array_get(h[code], i)]++;

  for (j = 1; j < ulen - q + 1; j++) {
    /* update code */
    code = qgram_step(code, u[j-1], u[j+q-1], alphabet_size,
                      r_raised_to_the_power_of_q_minus_1);
    /* increase count */
    for (i = 0; i < array_size(h[code]); i++)
      count[wlen - 1 + j - *(unsigned long*) array_get(h[code], i)]++;
  }

  /* determine score fasta */
  for (i = 0; i < ulen + wlen + 1; i++) {
    if (count[i] > maxcount)
      maxcount = count[i];
  }

  /* free space */
  ma_free(count);
  for (i = 0; i < hsize; i++)
    array_delete(h[i]);
  ma_free(h);

  return maxcount;
}
