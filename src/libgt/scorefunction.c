/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "scorefunction.h"
#include "xansi.h"

struct ScoreFunction {
  ScoreMatrix *sm;
  int deletion_score,
      insertion_score;
};

ScoreFunction* scorefunction_new(ScoreMatrix *sm, int deletion_score,
                                 int insertion_score, Env *env)
{
  ScoreFunction *sf;
  assert(sm);
  sf = env_ma_malloc(env, sizeof (ScoreFunction));
  sf->sm = sm;
  sf->deletion_score = deletion_score;
  sf->insertion_score = insertion_score;
  return sf;
}

int scorefunction_get_score(const ScoreFunction *s,
                            unsigned char idx1, unsigned char idx2)
{
  assert(s);
  return scorematrix_get_score(s->sm, idx1, idx2);
}

int scorefunction_get_deletion_score(const ScoreFunction *s)
{
  assert(s);
  return s->deletion_score;
}

int scorefunction_get_insertion_score(const ScoreFunction *s)
{
  assert(s);
  return s->insertion_score;
}

void scorefunction_delete(ScoreFunction *sf, Env *env)
{
  if (!sf) return;
  scorematrix_delete(sf->sm, env);
  env_ma_free(sf, env);
}
