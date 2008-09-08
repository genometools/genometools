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
#include "core/array.h"
#include "core/array2dim.h"
#include "core/parseutils.h"
#include "core/score_matrix.h"
#include "core/str.h"
#include "core/tokenizer.h"
#include "core/undef.h"
#include "core/xansi.h"

struct GT_ScoreMatrix {
  GT_Alpha *alpha;
  unsigned int dimension;
  int **scores;
};

GT_ScoreMatrix* gt_score_matrix_new(GT_Alpha *alpha)
{
  GT_ScoreMatrix *sm;
  assert(alpha);
  sm = gt_malloc(sizeof (GT_ScoreMatrix));
  sm->alpha = gt_alpha_ref(alpha);
  sm->dimension = gt_alpha_size(alpha);
  gt_array2dim_calloc(sm->scores, sm->dimension, sm->dimension);
  return sm;
}

static int parse_alphabet_line(GT_Array *index_to_gt_alpha_char_mapping,
                               Tokenizer *tz, GT_Error *err)
{
  GT_Str *token;
  char *tokenstr, amino_acid, parsed_characters[UCHAR_MAX] = { 0 };
  int had_err = 0;
  gt_error_check(err);
  assert(index_to_gt_alpha_char_mapping && tz);
  assert(!gt_array_size(index_to_gt_alpha_char_mapping));
  while ((token = tokenizer_get_token(tz))) {
    if (gt_str_length(token) > 2) {
      gt_error_set(err, "illegal character token '%s' on line %lu in file '%s'",
                gt_str_get(token), tokenizer_get_line_number(tz),
                tokenizer_get_filename(tz));
      had_err = -1;
      break;
    }
    tokenstr = gt_str_get(token);
    amino_acid = tokenstr[0];
    /* check for character duplications */
    if (parsed_characters[(int) amino_acid]) {
      gt_error_set(err, "the character '%c' appears more then once on line %lu in "
                "file  '%s'", amino_acid, tokenizer_get_line_number(tz),
                tokenizer_get_filename(tz));
      had_err = -1;
      break;
    }
    parsed_characters[(int) amino_acid] = UNDEF_CHAR;
    if (amino_acid == '\n') {
      gt_str_delete(token);
      tokenizer_next_token(tz);
      assert(!had_err);
      return 0;
    }
    gt_array_add(index_to_gt_alpha_char_mapping, amino_acid);
    if (gt_str_length(token) == 2) {
      if (tokenstr[1] != '\n') {
        gt_error_set(err, "illegal character token '%s' on line %lu in file '%s'",
                  gt_str_get(token), tokenizer_get_line_number(tz),
                  tokenizer_get_filename(tz));
        had_err = -1;
        break;
      }
      gt_str_delete(token);
      tokenizer_next_token(tz);
      assert(!had_err);
      return 0;
    }
    gt_str_delete(token);
    tokenizer_next_token(tz);
  }
  if (!had_err) {
    if (!gt_array_size(index_to_gt_alpha_char_mapping)) {
      gt_error_set(err, "could not parse a single alphabet character in file "
                "'%s' (file empty or directory?)", tokenizer_get_filename(tz));
    had_err = -1;
    }
  }
  gt_str_delete(token);
  return had_err;
}

static int parse_score_line(GT_ScoreMatrix *sm, Tokenizer *tz,
                            GT_Array *index_to_gt_alpha_char_mapping,
                            char *parsed_characters, GT_Error *err)
{
  unsigned int i = 0;
  char amino_acid;
  int score, had_err = 0;
  GT_Str *token;
  assert(sm && tz && index_to_gt_alpha_char_mapping);
  gt_error_check(err);
  token = tokenizer_get_token(tz);
  assert(token);
  if (gt_str_length(token) != 1) {
    gt_error_set(err, "illegal character token '%s' on line %lu in file '%s'",
              gt_str_get(token), tokenizer_get_line_number(tz),
              tokenizer_get_filename(tz));
    had_err = -1;
  }
  amino_acid = gt_str_get(token)[0];
  /* check for character duplications */
  if (parsed_characters[(int) amino_acid]) {
    gt_error_set(err, "multiple character '%c' entry on line %lu in file '%s'",
              amino_acid, tokenizer_get_line_number(tz),
              tokenizer_get_filename(tz));
    had_err = -1;
  }
  parsed_characters[(int) amino_acid] = UNDEF_CHAR;
  gt_str_delete(token);
  if (!had_err) {
    tokenizer_next_token(tz);
    while ((token = tokenizer_get_token(tz))) {
      had_err = gt_parse_int_line(&score, gt_str_get(token),
                                  tokenizer_get_line_number(tz),
                                  tokenizer_get_filename(tz), err);
      if (had_err)
        break;
      gt_score_matrix_set_score(sm,
                             gt_alpha_encode(sm->alpha, amino_acid),
                             gt_alpha_encode(sm->alpha, *(char*)
                             gt_array_get(index_to_gt_alpha_char_mapping, i)), score);
      i++;
      gt_str_delete(token);
      tokenizer_next_token(tz);
      if (tokenizer_line_start(tz))
          break;
    }
  }
  return had_err;
}

/* the score matrix parser */
static int parse_score_matrix(GT_ScoreMatrix *sm, const char *path, GT_Error *err)
{
  Tokenizer *tz;
  GT_Array *index_to_gt_alpha_char_mapping;
  unsigned int parsed_score_lines = 0;
  char parsed_characters[UCHAR_MAX] = { 0 };
  int had_err = 0;
  gt_error_check(err);
  assert(sm && path && sm->alpha);
  tz = tokenizer_new(gt_io_new(path, "r"));
  index_to_gt_alpha_char_mapping = gt_array_new(sizeof (char));
  tokenizer_skip_comment_lines(tz);
  had_err = parse_alphabet_line(index_to_gt_alpha_char_mapping, tz, err);
  if (!had_err) {
    while (tokenizer_has_token(tz)) {
      had_err = parse_score_line(sm, tz, index_to_gt_alpha_char_mapping,
                                 parsed_characters, err);
      if (had_err)
        break;
      parsed_score_lines++;
    }
  }

  /* check the number of parsed score lines */
  if (!had_err &&
      parsed_score_lines != gt_array_size(index_to_gt_alpha_char_mapping)) {
    gt_error_set(err, "the score matrix given in '%s' is not symmetric", path);
    had_err = -1;
  }

  gt_array_delete(index_to_gt_alpha_char_mapping);
  tokenizer_delete(tz);

  return had_err;
}

GT_ScoreMatrix* gt_score_matrix_new_read_protein(const char *path, GT_Error *err)
{
  GT_Alpha *protein_alpha;
  GT_ScoreMatrix *sm;
  int had_err;

  gt_error_check(err);
  assert(path);

  /* create score matrix */
  protein_alpha = gt_alpha_new_protein();
  sm = gt_score_matrix_new(protein_alpha);
  gt_alpha_delete(protein_alpha);

  /* parse matrix file */
  had_err = parse_score_matrix(sm, path, err);

  if (had_err) {
    gt_score_matrix_delete(sm);
    return NULL;
  }
  return sm;
}

unsigned int gt_score_matrix_get_dimension(const GT_ScoreMatrix *sm)
{
  assert(sm);
  return sm->dimension;
}

int gt_score_matrix_get_score(const GT_ScoreMatrix *sm,
                          unsigned int idx1, unsigned int idx2)
{
  assert(sm);
  assert(idx1 < sm->dimension && idx2 < sm->dimension); /* indices are valid */
  return sm->scores[idx1][idx2];
}

void gt_score_matrix_set_score(GT_ScoreMatrix *sm,
                           unsigned int idx1, unsigned int idx2, int score)
{
  assert(sm);
  assert(idx1 < sm->dimension && idx2 < sm->dimension); /* indices are valid */
  sm->scores[idx1][idx2] = score;
}

const int** gt_score_matrix_get_scores(const GT_ScoreMatrix *sm)
{
  assert(sm);
  return (const int**) sm->scores;
}

void gt_score_matrix_show(const GT_ScoreMatrix *sm, FILE *fp)
{
  unsigned i, j;
  assert(sm && fp);
  /* show alphabet line */
  xfputc(' ', fp);
  for (i = 0; i < gt_alpha_size(sm->alpha); i++)
    fprintf(fp, "  %c", gt_alpha_decode(sm->alpha, i));
  xfputc('\n', fp);
  /* show score lines */
  for (i = 0; i < gt_alpha_size(sm->alpha); i++) {
    xfputc(gt_alpha_decode(sm->alpha, i), fp);
    for (j = 0; j < gt_alpha_size(sm->alpha); j++)
      fprintf(fp, " %2d", gt_score_matrix_get_score(sm, i, j));
    xfputc('\n', fp);
  }
}

void gt_score_matrix_delete(GT_ScoreMatrix *sm)
{
  if (!sm) return;
  gt_alpha_delete(sm->alpha);
  gt_array2dim_delete(sm->scores);
  gt_free(sm);
}
