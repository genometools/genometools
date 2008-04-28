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
#include <limits.h>
#include "libgtcore/array2dim.h"
#include "libgtcore/minmax.h"
#include "libgtcore/undef.h"
#include "libgtext/coordinate.h"
#include "libgtext/swalign.h"

typedef struct {
  long score;
  bool max_replacement,
       max_deletion,
       max_insertion;
} DPentry;

static void fillDPtable(DPentry **dptable,
                        const char *u, unsigned long ulen,
                        const char *v, unsigned long vlen,
                        const int **scores,
                        int deletion_score, int insertion_score,
                        Coordinate *max_coordinate)
{
  unsigned long i, j;
  long maxscore, repscore, delscore, insscore, overall_maxscore = LONG_MIN;
  assert(dptable && u && ulen && v && vlen && max_coordinate);
  for (j = 1; j <= vlen; j++) {
    for (i = 1; i <= ulen; i++) {
      repscore = dptable[i-1][j-1].score + scores[(int) u[i-1]][(int) v[j-1]];
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

static Coordinate traceback(Alignment *a, DPentry **dptable,
                            unsigned long i, unsigned long j)
{
  Coordinate start_coordinate = { UNDEF_ULONG, UNDEF_ULONG };
  assert(a && dptable);
  while (dptable[i][j].score) {
    assert(dptable[i][j].score > 0);
    start_coordinate.x = i;
    start_coordinate.y = j;
    if (dptable[i][j].max_replacement) {
      alignment_add_replacement(a);
      i--;
      j--;
    }
    else if (dptable[i][j].max_deletion) {
      alignment_add_deletion(a);
      i--;;
    }
    else if (dptable[i][j].max_insertion) {
      alignment_add_insertion(a);
      j--;
    }
  }
  assert(start_coordinate.x != UNDEF_ULONG);
  assert(start_coordinate.y != UNDEF_ULONG);
  return start_coordinate;
}

static Alignment* smith_waterman_align(const char *u_orig, const char *v_orig,
                                       const char *u_enc, const char *v_enc,
                                       unsigned long u_len, unsigned long v_len,
                                       const int **scores,
                                       int deletion_score, int insertion_score)
{
  assert(u_orig && v_orig && u_enc && v_enc && u_len && v_len && scores);
  Coordinate alignment_start, alignment_end = { UNDEF_ULONG, UNDEF_ULONG };
  DPentry **dptable;
  Alignment *a = NULL;
  array2dim_calloc(dptable, u_len+1, v_len+1);
  fillDPtable(dptable, u_enc, u_len, v_enc, v_len, scores, deletion_score,
              insertion_score, &alignment_end);
  assert(alignment_end.x != UNDEF_ULONG);
  assert(alignment_end.y != UNDEF_ULONG);
  if (dptable[alignment_end.x][alignment_end.y].score) {
    /* construct only an alignment if a (positive) score was computed */
    a = alignment_new();
    alignment_start = traceback(a, dptable, alignment_end.x, alignment_end.y);
    /* transform the positions in the DP matrix to sequence positions */
    alignment_start.x--;
    alignment_start.y--;
    alignment_end.x--;
    alignment_end.y--;
    /* employ sequence positions to set alignment sequences */
    alignment_set_seqs(a,
                       u_orig + alignment_start.x,
                       alignment_end.x - alignment_start.x + 1,
                       v_orig + alignment_start.y,
                       alignment_end.y - alignment_start.y + 1);
  }
  array2dim_delete(dptable);
  return a;
}

Alignment* swalign(Seq *u, Seq *v, const ScoreFunction *sf)
{
  assert(u && v && sf);
  return smith_waterman_align(seq_get_orig(u), seq_get_orig(v),
                              seq_get_encoded(u), seq_get_encoded(v),
                              seq_length(u), seq_length(v),
                              score_function_get_scores(sf),
                              score_function_get_deletion_score(sf),
                              score_function_get_insertion_score(sf));
}
