/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgt/gc_content.h>

void gc_content_show(const char *seq, unsigned long len, Alpha *alpha, Env *env)
{
  unsigned long i,
                gc = 0, /* number of G/C bases */
                at = 0, /* number of A/T bases */
                n  = 0; /* number of N   bases */
  unsigned int a_code, c_code, g_code, t_code, n_code, cc;
  Alpha *dna_alpha;
  env_error_check(env);
  assert(seq && alpha);
  dna_alpha = alpha_new_dna(env);
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
  alpha_delete(dna_alpha, env);
}
