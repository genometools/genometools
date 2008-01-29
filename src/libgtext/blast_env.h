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

#ifndef BLAST_ENV_H
#define BLAST_ENV_H

#include "libgtcore/score_matrix.h"

typedef struct BlastEnv BlastEnv;

/*
   Construct the BlastP environment for encoded sequence <w> of length <wlen>
   (over alphabet <alpha>).  Thereby, <q> is the q-gram length and <k> is the
   minimum score. The <score_matrix> is used to determine the match/mismatch
   scores.
   Returns a new BlastEnv object.
*/
BlastEnv* blast_env_new(const char *w, unsigned long wlen, Alpha *alpha,
                        unsigned long q, long k,
                        const ScoreMatrix *score_matrix);

/* Delete the Blast environment <blast_env>. */
void      blast_env_delete(BlastEnv*);

/* Show the Blast environment <blast_env> on stdout. */
void      blast_env_show(const BlastEnv*);

#endif
