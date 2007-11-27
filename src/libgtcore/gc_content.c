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

#include <assert.h>
#include "libgtcore/gc_content.h"

void gc_content_show(const char *seq, unsigned long len, Alpha *alpha)
{
  unsigned long i,
                gc = 0, /* number of G/C bases */
                at = 0, /* number of A/T bases */
                n  = 0; /* number of N   bases */
  unsigned int a_code, c_code, g_code, t_code, n_code, cc;
  Alpha *dna_alpha;
  assert(seq && alpha);
  dna_alpha = alpha_new_dna();
  assert(alpha_is_compatible_with_alpha(alpha, dna_alpha));
  a_code = alpha_encode(dna_alpha, 'A');
  c_code = alpha_encode(dna_alpha, 'C');
  g_code = alpha_encode(dna_alpha, 'G');
  t_code = alpha_encode(dna_alpha, 'T');
  n_code = alpha_encode(dna_alpha, 'N');
  for (i = 0; i < len; i++) {
    cc = alpha_encode(alpha, seq[i]);
    if (cc == g_code || cc == c_code)
      gc++;
    else if (cc == a_code || cc == t_code)
      at++;
    else if (cc == n_code)
      n++;
    else {
      assert(0);
    }
  }
  printf("GC-content: %.2f%% (AT-content: %.2f%%, N-content: %.2f%%)\n",
         ((double) gc / len) * 100.0, ((double) at / len) * 100.0,
         ((double) n  / len) * 100.0);
  alpha_delete(dna_alpha);
}
