/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef SCORE_MATRIX_H
#define SCORE_MATRIX_H

#include "core/alphabet.h"
#include "core/error.h"

typedef struct GtScoreMatrix GtScoreMatrix;

/* A score matrix is always defined over a given <alphabet>. */
GtScoreMatrix* gt_score_matrix_new(GtAlphabet *alphabet);
/* Read in a protein score matrix from the given <path> and return it. */
GtScoreMatrix* gt_score_matrix_new_read_protein(const char *path, GtError *err);
/* Read in score matrix from <path> over given <alphabet> and return it. */
GtScoreMatrix* gt_score_matrix_new_read(const char *path, GtAlphabet *alphabet,
                                        GtError *err);
unsigned int   gt_score_matrix_get_dimension(const GtScoreMatrix*);
int            gt_score_matrix_get_score(const GtScoreMatrix*,
                                         unsigned int, unsigned int);
void           gt_score_matrix_set_score(GtScoreMatrix*,
                                         unsigned int, unsigned int, int);
const int**    gt_score_matrix_get_scores(const GtScoreMatrix*);
void           gt_score_matrix_show(const GtScoreMatrix*, FILE*);
void           gt_score_matrix_delete(GtScoreMatrix*);

#endif
