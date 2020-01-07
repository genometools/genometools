/*
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
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

#ifndef SCORE_MATRIX_API_H
#define SCORE_MATRIX_API_H

#include "core/alphabet_api.h"
#include "core/error_api.h"

/* <GtScoreMatrix> represents a matrix of signed integer values, for example
   for scoring alphabet symbols against each other. */
typedef struct GtScoreMatrix GtScoreMatrix;

/* A score matrix is always defined over a given <alphabet>. */
GtScoreMatrix* gt_score_matrix_new(GtAlphabet *alphabet);
/* Create empty score matrix with same dimension as <scorematrix>. */
GtScoreMatrix* gt_score_matrix_clone_empty(const GtScoreMatrix *scorematrix);
/* Read in a protein score matrix from the given <path> and return it. */
GtScoreMatrix* gt_score_matrix_new_read_protein(const char *path, GtError *err);
/* Read in score matrix from <path> over given <alphabet> and return it. */
GtScoreMatrix* gt_score_matrix_new_read(const char *path, GtAlphabet *alphabet,
                                        GtError *err);
/* Return the dimension of <scorematrix>. */
unsigned int   gt_score_matrix_get_dimension(const GtScoreMatrix *scorematrix);
/* Return the score value in <scorematrix> at positions (<idx1>,<idx2>). */
int            gt_score_matrix_get_score(const GtScoreMatrix *scorematrix,
                                         unsigned int idx1, unsigned int idx2);
/* Set the score value in <scorematrix> at positions (<idx1>,<idx2>) to
   <score>. */
void           gt_score_matrix_set_score(GtScoreMatrix *scorematrix,
                                         unsigned int idx1, unsigned int idx2,
                                         int score);
/* Return the score values in <scorematrix> as a two-dimensional array. */
const int**    gt_score_matrix_get_scores(const GtScoreMatrix *scorematrix);
/* Print <scorematrix> to output <fp>. */
void           gt_score_matrix_show(const GtScoreMatrix *scorematrix, FILE *fp);
/* Delete <scorematrix>. */
void           gt_score_matrix_delete(GtScoreMatrix *scorematrix);

#endif
