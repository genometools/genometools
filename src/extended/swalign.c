/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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
#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/minmax.h"
#include "core/undef_api.h"
#include "extended/swalign.h"

typedef struct {
  unsigned long x,
                y;
} Coordinate;

typedef struct {
  long score;
  bool max_replacement,
       max_deletion,
       max_insertion;
} DPentry;

static void swalign_fill_table(DPentry **dptable,
                               const GtUchar *u, unsigned long ulen,
                               const GtUchar *v, unsigned long vlen,
                               const int **scores,
                               int deletion_score, int insertion_score,
                               Coordinate *max_coordinate,
                               unsigned int u_alpha_size,
                               unsigned int v_alpha_size)
{
  unsigned long i, j;
  long maxscore, repscore, delscore, insscore, overall_maxscore = LONG_MIN;
  gt_assert(dptable && u && ulen && v && vlen && max_coordinate && u_alpha_size
            && v_alpha_size);
  for (j = 1; j <= vlen; j++) {
    for (i = 1; i <= ulen; i++) {
      int uval, vval;
      uval = (int) ((u[i-1] == WILDCARD) ? u_alpha_size - 1 : u[i-1]);
      vval = (int) ((v[j-1] == WILDCARD) ? v_alpha_size - 1 : v[j-1]);
      repscore = dptable[i-1][j-1].score + scores[uval][vval];
      delscore = dptable[i-1][j].score + deletion_score;
      insscore = dptable[i][j-1].score + insertion_score;
      maxscore = MAX(MAX(MAX(repscore, delscore), insscore), 0);
      dptable[i][j].score = maxscore;
      dptable[i][j].max_replacement = (maxscore == repscore) ? true : false;
      dptable[i][j].max_deletion    = (maxscore == delscore) ? true : false;
      dptable[i][j].max_insertion   = (maxscore == insscore) ? true : false;
      if (maxscore > overall_maxscore) {
        overall_maxscore = maxscore;
        max_coordinate->x = i;
        max_coordinate->y = j;
      }
    }
  }
}

static Coordinate traceback(GtAlignment *a, DPentry **dptable,
                            unsigned long i, unsigned long j)
{
  Coordinate start_coordinate = { GT_UNDEF_ULONG, GT_UNDEF_ULONG };
  gt_assert(a && dptable);
  while (dptable[i][j].score) {
    gt_assert(dptable[i][j].score > 0);
    start_coordinate.x = i;
    start_coordinate.y = j;
    if (dptable[i][j].max_replacement) {
      gt_alignment_add_replacement(a);
      i--;
      j--;
    }
    else if (dptable[i][j].max_deletion) {
      gt_alignment_add_deletion(a);
      i--;;
    }
    else if (dptable[i][j].max_insertion) {
      gt_alignment_add_insertion(a);
      j--;
    }
  }
  gt_assert(start_coordinate.x != GT_UNDEF_ULONG);
  gt_assert(start_coordinate.y != GT_UNDEF_ULONG);
  return start_coordinate;
}

static GtAlignment* smith_waterman_align(const char *u_orig,
                                         const char *v_orig,
                                         const GtUchar *u_enc,
                                         const GtUchar *v_enc,
                                         unsigned long u_len,
                                         unsigned long v_len,
                                         const int **scores,
                                         int deletion_score,
                                         int insertion_score,
                                         const GtAlphabet *u_alpha,
                                         const GtAlphabet *v_alpha)
{
  gt_assert(u_orig && v_orig && u_enc && v_enc && u_len && v_len && scores
            && u_alpha && v_alpha);
  Coordinate alignment_start,
             alignment_end = { GT_UNDEF_ULONG, GT_UNDEF_ULONG };
  GtRange urange, vrange;
  DPentry **dptable;
  GtAlignment *a = NULL;
  gt_array2dim_calloc(dptable, u_len+1, v_len+1);
  swalign_fill_table(dptable, u_enc, u_len, v_enc, v_len, scores,
                     deletion_score, insertion_score, &alignment_end,
                     gt_alphabet_size(u_alpha), gt_alphabet_size(v_alpha));
  gt_assert(alignment_end.x != GT_UNDEF_ULONG);
  gt_assert(alignment_end.y != GT_UNDEF_ULONG);
  if (dptable[alignment_end.x][alignment_end.y].score) {
    /* construct only an alignment if a (positive) score was computed */
    a = gt_alignment_new();
    alignment_start = traceback(a, dptable, alignment_end.x, alignment_end.y);
    /* transform the positions in the DP matrix to sequence positions */
    urange.start = --alignment_start.x;
    vrange.start = --alignment_start.y;
    urange.end = --alignment_end.x;
    vrange.end = --alignment_end.y;
    /* employ sequence positions to set alignment sequences */
    gt_alignment_set_seqs(a,
                          (const GtUchar *) (u_orig + alignment_start.x),
                          alignment_end.x - alignment_start.x + 1,
                          (const GtUchar *) (v_orig + alignment_start.y),
                          alignment_end.y - alignment_start.y + 1);
    gt_alignment_set_urange(a, urange);
    gt_alignment_set_vrange(a, vrange);
  }
  gt_array2dim_delete(dptable);
  return a;
}

GtAlignment* gt_swalign(GtSeq *u, GtSeq *v, const GtScoreFunction *sf)
{
  gt_assert(u && v && sf);
  return smith_waterman_align(gt_seq_get_orig(u), gt_seq_get_orig(v),
                              gt_seq_get_encoded(u), gt_seq_get_encoded(v),
                              gt_seq_length(u), gt_seq_length(v),
                              gt_score_function_get_scores(sf),
                              gt_score_function_get_deletion_score(sf),
                              gt_score_function_get_insertion_score(sf),
                              gt_seq_get_alphabet(u), gt_seq_get_alphabet(v));
}
