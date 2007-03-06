/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <limits.h>
#include <libgt/array2dim.h>
#include <libgt/coordinate.h>
#include <libgt/minmax.h>
#include <libgt/swalign.h>

typedef struct {
  long score;
  bool max_replacement,
       max_deletion,
       max_insertion;
} DPentry;

static void fillDPtable(DPentry **dptable,
                        const char *u, unsigned long ulen,
                        const char *v, unsigned long vlen,
                        const ScoreFunction *sf,
                        Coordinate *max_coordinate)
{
  unsigned long i, j;
  long maxscore, repscore, delscore, insscore, overall_maxscore = LONG_MIN;
  assert(dptable && u && ulen && v && vlen && max_coordinate);
  for (j = 1; j <= vlen; j++) {
    for (i = 1; i <= ulen; i++) {
      repscore = dptable[i-1][j-1].score +
                 scorefunction_get_score(sf, u[i-1], v[j-1]);
      delscore = dptable[i-1][j].score + scorefunction_get_deletion_score(sf);
      insscore = dptable[i][j-1].score + scorefunction_get_insertion_score(sf);
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
                            unsigned long i, unsigned long j, Env *env)
{
  Coordinate start_coordinate;
  assert(a && dptable);
  while (dptable[i][j].score) {
    assert(dptable[i][j].score > 0);
    start_coordinate.x = i;
    start_coordinate.y = j;
    if (dptable[i][j].max_replacement) {
      alignment_add_replacement(a, env);
      i--;
      j--;
    }
    else if (dptable[i][j].max_deletion) {
      alignment_add_deletion(a, env);
      i--;;
    }
    else if (dptable[i][j].max_insertion) {
      alignment_add_insertion(a, env);
      j--;
    }
  }
  return start_coordinate;
}

Alignment* swalign(Seq *u, Seq *v, const ScoreFunction *sf, Env *env)
{
  Coordinate alignment_start, alignment_end;
  DPentry **dptable;
  Alignment *a = NULL;
  assert(u && v && sf);
  array2dim_calloc(dptable, seq_length(u)+1, seq_length(v)+1, DPentry, env);
  fillDPtable(dptable, seq_get_encoded(u, env), seq_length(u),
              seq_get_encoded(v, env), seq_length(v), sf, &alignment_end);
  if (dptable[alignment_end.x][alignment_end.y].score) {
    /* construct only an alignment if a (positive) score was computed */
    a = alignment_new(env);
    alignment_start = traceback(a, dptable, alignment_end.x, alignment_end.y,
                                env);
    /* transform the positions in the DP matrix to sequence positions */
    alignment_start.x--;
    alignment_start.y--;
    alignment_end.x--;
    alignment_end.y--;
    /* employ sequence positions to set alignment sequences */
    alignment_set_seqs(a,
                       seq_get_orig(u) + alignment_start.x,
                       alignment_end.x - alignment_start.x + 1,
                       seq_get_orig(v) + alignment_start.y,
                       alignment_end.y - alignment_start.y + 1);
  }
  array2dim_delete(dptable, env);
  return a;
}
