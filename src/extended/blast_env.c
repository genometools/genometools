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
#include "core/bittab.h"
#include "core/hashmap-generic.h"
#include "core/ma.h"
#include "core/xansi.h"
#include "extended/blast_env.h"
#include "extended/qgram.h"

/* Internal class to hold position lists. */
typedef struct {
  GtHashtable *mapping;
} GtPos;

DECLARE_HASHMAP(unsigned long, ul, GtArray *, array, static, inline)
DEFINE_HASHMAP(unsigned long, ul, GtArray *, array, gt_ht_ul_elem_hash,
               gt_ht_ul_elem_cmp, NULL_DESTRUCTOR, gt_array_delete, static,
               inline)
DECLARE_SAFE_DEREF(GtArray *,array)

/* Return a new GtPos object. */
GtPos* gt_pos_new(void)
{
  GtPos *pos = gt_malloc(sizeof *pos);
  pos->mapping = ul_array_gt_hashmap_new();
  return pos;
}

/* Delete <pos> object. */
void gt_pos_delete(GtPos *pos)
{
  if (!pos) return;
  gt_hashtable_delete(pos->mapping);
  gt_free(pos);
}

/* Add <position> to position list for <code>. */
void gt_pos_add(GtPos *pos, unsigned long code, unsigned long position)
{
  GtArray *position_list;
  gt_assert(pos && pos->mapping);
  position_list = array_gt_safe_deref(ul_array_gt_hashmap_get(pos->mapping,
                                      code));
  if (!position_list) {
    position_list = gt_array_new(sizeof (unsigned long));
    gt_array_add(position_list, position);
    ul_array_gt_hashmap_add(pos->mapping, code, position_list);
  }
  else
    gt_array_add(position_list, position);
}

/* Get position list for <code>. */
GtArray* gt_pos_get(GtPos *pos, unsigned long code)
{
  gt_assert(pos && pos->mapping);
  return array_gt_safe_deref(ul_array_gt_hashmap_get(pos->mapping, code));
}

/*
  Compute the maximal position scores for sequence <w> of length <wlen> under
  score matrix <score_matrix>. A maximal position score is the maximal score
  which a certain position of <w> can contribute to the total score.
  Returns an array of length <wlen> which has to be freed by the caller.
*/
static long* compute_max_pos_scores(const char *w, unsigned long wlen,
                                    const GtScoreMatrix *score_matrix)
{
  long score, max_score, *max_matrix_scores, *max_pos_scores;
  unsigned int dimension;
  unsigned long i;

  dimension = gt_score_matrix_get_dimension(score_matrix);
  max_matrix_scores = gt_malloc(sizeof (long) * dimension);
  max_pos_scores = gt_malloc(sizeof (long) * wlen);
  /* fill maximal matrix scores */
  for (i = 0; i < dimension; i++) {
    unsigned long j;
    max_score = LONG_MIN;
    for (j = 0; j < dimension; j++) {
      score = gt_score_matrix_get_score(score_matrix, i, j);
      if (score > max_score)
        max_score = score;
    }
    max_matrix_scores[i] = max_score;
  }
  /* fill maximal position scores */
  for (i = 0; i < wlen; i++)
    max_pos_scores[i] = max_matrix_scores[(int) w[i]];
  /* free */
  gt_free(max_matrix_scores);

  return max_pos_scores;
}

/*
  This function does the actual work of adding all strings with score >= <k> to
  the Blast environment (which consists of the bit vector <V> and the set of
  position lists <pos>).
  This is achieved by constructing an implicit trie (via recursive
  function calls) representing all possible strings of length <q> over <alpha>.
  - <qgram_rest> is the rest of the current q-gram for which the trie is
    constructed.
  - <current_word> holds the word which denotes the characters along
    the path to the current node.
  - <max_cumul_scores> contains a the cumulative scores necessary for the
    lookahead mechanism.
  - <alpha> is the alphabet.
  - <q> is the q-gram size.
  - <q_rest> is the size of <qgram_rest>.
  - <k> is the minimum score which is necessary to add the <current_word> to the
    Blast environment).
  - <current_score> denotes the score of the <current_word>.
  - <position> denotes the position of the current q-gram in the query string.

  The lookahead mechanism works as follows:
  If the <current_score> of the <current_word> plus the maximal score achievable
  during the rest of the traversal (stored in <max_cumul_scores>) is smaller
  than the minimum score <k>, then the traversal of this branch is aborted
  (because it is not possible to find a word with a sufficient score on this
  branch).
*/
static void add_q_word_to_env(GtBittab *V, GtPos *pos, const char *qgram_rest,
                              char *current_word, long *max_cumul_scores,
                              GtAlpha *alpha, unsigned long q,
                              unsigned long q_rest, long k, long current_score,
                              unsigned long position,
                              const GtScoreMatrix *score_matrix)
{
  gt_assert(V && pos && qgram_rest && alpha);
  if (q_rest == 0) {
    if (current_score >= k) {
      unsigned long qgram_code = gt_qgram_encode(current_word, q,
                                                 gt_alpha_size(alpha));
      /* set V[qgram_code] */
      gt_bittab_set_bit(V, qgram_code);
      /* store position */
      gt_pos_add(pos, qgram_code, position);
    }
  }
  else if (current_score + max_cumul_scores[q - q_rest] >= k) { /* lookahead */
    unsigned int i;
    /* enumerate all possible characters (at this node) */
    for (i = 0; i < gt_alpha_size(alpha); i++) {
      long char_score = gt_score_matrix_get_score(score_matrix, qgram_rest[0],
                                                  i);
      gt_assert(i < CHAR_MAX);
      current_word[q - q_rest] = i;
      /* recursive call (descend into branch) */
      add_q_word_to_env(V, pos, qgram_rest+1, current_word, max_cumul_scores,
                        alpha, q, q_rest-1, k, current_score+char_score,
                        position, score_matrix);
    }
  }
}

/*
  This function fills the Blast environment, which consists of the bit vector
  <V> and the set of position lists <pos>, for the encoded word <w> of length
  <wlen>.
*/
static void compute_env(GtBittab *V, GtPos *pos, const char *w,
                        unsigned long wlen, GtAlpha *alpha, unsigned long q,
                        long k, const GtScoreMatrix *score_matrix)
{
  long *max_pos_scores, *max_cumul_scores;
  unsigned long i;
  char *current_word;
  /* prepare space for current word */
  current_word = gt_malloc(sizeof (char) * (q+1));
  current_word[q] = '\0';
  /* prepare data for lookahead */
  max_pos_scores = compute_max_pos_scores(w, wlen, score_matrix);
  max_cumul_scores = gt_malloc(sizeof (long) * q);
  /* add all words to <blast_env> */
  for (i = 0; i < wlen - q + 1; i++) {
    long j;
    /* fill maximal cumulative scores (for lookahead) */
    max_cumul_scores[q-1] = max_pos_scores[i+q-1];
    for (j = q-2; j >= 0; j--)
      max_cumul_scores[j] = max_pos_scores[i+j] + max_cumul_scores[j+1];
    /* process the current q-word */
    add_q_word_to_env(V, pos, w+i, current_word, max_cumul_scores, alpha, q, q,
                      k, 0, i, score_matrix);
  }
  /* free */
  gt_free(max_cumul_scores);
  gt_free(max_pos_scores);
  gt_free(current_word);
}

/* the actual class implementation */
struct GtBlastEnv {
  GtAlpha *alpha;
  unsigned long q;
  GtBittab *V; /* The vector V of r^q bits. */
  GtPos *pos;  /* The set of position lists. If a bit in <V> is set then <pos>
                contains the corresponding position list for that code. */
};

GtBlastEnv* gt_blast_env_new(const char *w, unsigned long wlen, GtAlpha *alpha,
                             unsigned long q, long k,
                             const GtScoreMatrix *score_matrix)
{
  GtBlastEnv *be;
  gt_assert(w && alpha && q && score_matrix);
  gt_assert(gt_alpha_size(alpha) ==
            gt_score_matrix_get_dimension(score_matrix));
  be = gt_calloc(1, sizeof *be);
  be->alpha = gt_alpha_ref(alpha);
  be->q = q;
  /* if <w> is long enough fill the Blast environment */
  if (wlen >= q) {
    be->V = gt_bittab_new(pow(gt_alpha_size(alpha), q));
    be->pos = gt_pos_new();
    compute_env(be->V, be->pos, w, wlen, alpha, q, k, score_matrix);
  }
  return be;
}

void gt_blast_env_delete(GtBlastEnv *be)
{
  if (!be) return;
  gt_alpha_delete(be->alpha);
  gt_pos_delete(be->pos);
  gt_bittab_delete(be->V);
  gt_free(be);
}

void gt_blast_env_show(const GtBlastEnv *be)
{
  unsigned long i, code;
  GtArray *position_list;
  char *qgram;
  gt_assert(be);
  if (!be->pos)
    return; /* nothing to do */
  gt_assert(be->V && be->alpha);
  qgram = gt_malloc(sizeof (char) * (be->q+1));
  qgram[be->q] = '\0';
  /* iterate over all codes which have a position list */
  for (code  = gt_bittab_get_first_bitnum(be->V);
       code != gt_bittab_get_last_bitnum(be->V);
       code  = gt_bittab_get_next_bitnum(be->V, code)) {
    position_list = gt_pos_get(be->pos, code);
    gt_assert(position_list);
    gt_assert(gt_array_size(position_list)); /* contains at least one pos. */
    gt_qgram_decode(qgram, code, be->q, be->alpha);
    gt_xfputs(qgram, stdout);
    for (i = 0; i < gt_array_size(position_list); i++) {
      printf(", %lu", *(unsigned long*) gt_array_get(position_list, i) + 1);
    }
    gt_xputchar('\n');
  }
  gt_free(qgram);
}
