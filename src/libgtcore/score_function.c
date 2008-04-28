/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/ma.h"
#include "libgtcore/score_function.h"
#include "libgtcore/xansi.h"

struct ScoreFunction {
  ScoreMatrix *sm;
  int deletion_score,
      insertion_score;
};

ScoreFunction* score_function_new(ScoreMatrix *sm, int deletion_score,
                                 int insertion_score)
{
  ScoreFunction *sf;
  assert(sm);
  sf = ma_malloc(sizeof (ScoreFunction));
  sf->sm = sm;
  sf->deletion_score = deletion_score;
  sf->insertion_score = insertion_score;
  return sf;
}

int score_function_get_score(const ScoreFunction *sf,
                            unsigned int idx1, unsigned int idx2)
{
  assert(sf);
  return score_matrix_get_score(sf->sm, idx1, idx2);
}

const int** score_function_get_scores(const ScoreFunction *sf)
{
  assert(sf);
  return score_matrix_get_scores(sf->sm);
}

int score_function_get_deletion_score(const ScoreFunction *sf)
{
  assert(sf);
  return sf->deletion_score;
}

int score_function_get_insertion_score(const ScoreFunction *sf)
{
  assert(sf);
  return sf->insertion_score;
}

void score_function_delete(ScoreFunction *sf)
{
  if (!sf) return;
  score_matrix_delete(sf->sm);
  ma_free(sf);
}
