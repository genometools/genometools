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

#include "core/ma.h"
#include "core/score_function.h"

struct GtScoreFunction {
  GtScoreMatrix *sm;
  int deletion_score,
      insertion_score;
  unsigned long reference_count;
};

GtScoreFunction* gt_score_function_new(GtScoreMatrix *sm, int deletion_score,
                                 int insertion_score)
{
  GtScoreFunction *sf;
  gt_assert(sm);
  sf = gt_malloc(sizeof (GtScoreFunction));
  sf->sm = sm;
  sf->deletion_score = deletion_score;
  sf->insertion_score = insertion_score;
  sf->reference_count = 0;
  return sf;
}

GtScoreFunction* gt_score_function_ref(GtScoreFunction *sf)
{
  if (!sf) return NULL;
  sf->reference_count++;
  return sf;
}

int gt_score_function_get_score(const GtScoreFunction *sf,
                            unsigned int idx1, unsigned int idx2)
{
  gt_assert(sf);
  return gt_score_matrix_get_score(sf->sm, idx1, idx2);
}

const int** gt_score_function_get_scores(const GtScoreFunction *sf)
{
  gt_assert(sf);
  return gt_score_matrix_get_scores(sf->sm);
}

int gt_score_function_get_deletion_score(const GtScoreFunction *sf)
{
  gt_assert(sf);
  return sf->deletion_score;
}

int gt_score_function_get_insertion_score(const GtScoreFunction *sf)
{
  gt_assert(sf);
  return sf->insertion_score;
}

void gt_score_function_delete(GtScoreFunction *sf)
{
  if (!sf) return;
  if (sf->reference_count) {
    sf->reference_count--;
    return;
  }
  gt_score_matrix_delete(sf->sm);
  gt_free(sf);
}
