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

#include <limits.h>
#include <math.h>
#include "libgtcore/bittab.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/xansi.h"
#include "libgtext/blast_env.h"
#include "libgtext/qgram.h"

typedef struct{
  Hashtable *mapping;
} Pos;

struct BlastEnv {
  Alpha *alpha;
  unsigned long q;
  Bittab *V; /* the vector V of r^q bits */
  Pos *pos; /* the set of position lists */
};

Pos* pos_new(void)
{
  Pos *pos = ma_malloc(sizeof *pos);
  pos->mapping = hashtable_new(HASH_DIRECT, NULL, (FreeFunc) array_delete);
  return pos;
}

void pos_delete(Pos *pos)
{
  if (!pos) return;
  hashtable_delete(pos->mapping);
  ma_free(pos);
}

void pos_add(Pos *pos, unsigned long code, unsigned long position)
{
  Array *position_list;
  assert(pos && pos->mapping);
  position_list = hashtable_get(pos->mapping, (void*) code);
  if (!position_list) {
    position_list = array_new(sizeof (unsigned long));
    array_add(position_list, position);
    hashtable_add(pos->mapping, (void*) code, position_list);
  }
  else
    array_add(position_list, position);
}

Array* pos_get(Pos *pos, unsigned long code)
{
  assert(pos && pos->mapping);
  return hashtable_get(pos->mapping, (void*) code);
}

static void blast_env_add_q_word(BlastEnv *be, const char *qgram_rest,
                                 char *current_word, Alpha *alpha,
                                 unsigned long q, unsigned long q_rest, long k,
                                 long current_score, unsigned long position,
                                 const ScoreMatrix *score_matrix)
{
  assert(be && qgram_rest && alpha);
  if (q_rest == 0) {
    if (current_score >= k) {
      unsigned long qgram_code = qgram_encode(current_word, q,
                                              alpha_size(alpha));
      /* set V[qgram_code] */
      bittab_set_bit(be->V, qgram_code);
      /* store position */
      pos_add(be->pos, qgram_code, position);
    }
  }
  else {
    unsigned int i;
    for (i = 0; i < alpha_size(alpha); i++) {
      long char_score = score_matrix_get_score(score_matrix, qgram_rest[0], i);
      assert(i < CHAR_MAX);
      current_word[q - q_rest] = i;
      blast_env_add_q_word(be, qgram_rest+1, current_word, alpha, q, q_rest-1,
                           k, current_score+char_score, position, score_matrix);
    }
  }
}

BlastEnv* blast_env_new(const char *w, unsigned long wlen, Alpha *alpha,
                        unsigned long q, long k,
                        const ScoreMatrix *score_matrix)
{
  unsigned int r; /* the alphabet size */
  BlastEnv *be;
  assert(w && alpha && q && score_matrix);
  r = alpha_size(alpha);
  assert(r == score_matrix_get_dimension(score_matrix));
  be = ma_calloc(1, sizeof *be);
  be->alpha = alpha_ref(alpha);
  be->q = q;
  if (wlen >= q) {
    unsigned long i;
    char *current_word;
    be->V = bittab_new(pow(r, q));
    be->pos = pos_new();
    /* prepare space for current word */
    current_word = ma_malloc(sizeof (char) * (q+1));
    current_word[q] = '\0';
    /* add all words to <blast_env> */
    for (i = 0; i < wlen - q + 1; i++) {
      blast_env_add_q_word(be, w+i, current_word, alpha, q, q, k, 0, i,
                           score_matrix);
    }
    /* free */
    ma_free(current_word);
  }
  return be;
}

void blast_env_delete(BlastEnv *be)
{
  if (!be) return;
  alpha_delete(be->alpha);
  pos_delete(be->pos);
  bittab_delete(be->V);
  ma_free(be);
}

void blast_env_show(const BlastEnv *be)
{
  unsigned long i, code;
  Array *position_list;
  char *qgram;
  assert(be);
  if (!be->pos)
    return; /* nothing to do */
  assert(be->V && be->alpha);
  qgram = ma_malloc(sizeof (char) * (be->q+1));
  qgram[be->q] = '\0';
  /* iterate over all codes which have a position list */
  for (code  = bittab_get_first_bitnum(be->V);
       code != bittab_get_last_bitnum(be->V);
       code  = bittab_get_next_bitnum(be->V, code)) {
    position_list = pos_get(be->pos, code);
    assert(position_list);
    assert(array_size(position_list)); /* contains at least one position */
    qgram_decode(qgram, code, be->q, be->alpha);
    xfputs(qgram, stdout);
    for (i = 0; i < array_size(position_list); i++) {
      printf(", %lu", *(unsigned long*) array_get(position_list, i) + 1);
    }
    xputchar('\n');
  }
  ma_free(qgram);
}
