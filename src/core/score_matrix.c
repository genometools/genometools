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

#include <string.h>
#include "core/assert_api.h"
#include "core/array.h"
#include "core/array2dim_api.h"
#include "core/chardef.h"
#include "core/parseutils.h"
#include "core/score_matrix.h"
#include "core/str.h"
#include "core/tokenizer.h"
#include "core/undef_api.h"
#include "core/xansi_api.h"

struct GtScoreMatrix {
  GtAlphabet *alphabet;
  unsigned int dimension;
  int **scores;
};

GtScoreMatrix* gt_score_matrix_new(GtAlphabet *alphabet)
{
  GtScoreMatrix *sm;
  gt_assert(alphabet);
  sm = gt_malloc(sizeof (GtScoreMatrix));
  sm->alphabet = gt_alphabet_ref(alphabet);
  sm->dimension = gt_alphabet_size(alphabet);
  gt_array2dim_calloc(sm->scores, sm->dimension, sm->dimension);
  return sm;
}

static int parse_alphabet_line(GtArray *index_to_alpha_char_mapping,
                               GtTokenizer *tz, GtError *err)
{
  GtStr *token;
  char *tokenstr, amino_acid, parsed_characters[UCHAR_MAX] = { 0 };
  int had_err = 0;
  gt_error_check(err);
  gt_assert(index_to_alpha_char_mapping && tz);
  gt_assert(!gt_array_size(index_to_alpha_char_mapping));
  while ((token = gt_tokenizer_get_token(tz))) {
    if (gt_str_length(token) > 2) {
      gt_error_set(err, "illegal character token '%s' on line %lu in file '%s'",
                gt_str_get(token), gt_tokenizer_get_line_number(tz),
                gt_tokenizer_get_filename(tz));
      had_err = -1;
      break;
    }
    tokenstr = gt_str_get(token);
    amino_acid = tokenstr[0];
    /* check for character duplications */
    if (parsed_characters[(int) amino_acid]) {
      gt_error_set(err, "the character '%c' appears more then once on line %lu "
                   "in file  '%s'", amino_acid,
                   gt_tokenizer_get_line_number(tz),
                   gt_tokenizer_get_filename(tz));
      had_err = -1;
      break;
    }
    parsed_characters[(int) amino_acid] = GT_UNDEF_CHAR;
    if (amino_acid == '\n') {
      gt_str_delete(token);
      gt_tokenizer_next_token(tz);
      gt_assert(!had_err);
      return 0;
    }
    gt_array_add(index_to_alpha_char_mapping, amino_acid);
    if (gt_str_length(token) == 2) {
      if (tokenstr[1] != '\n') {
        gt_error_set(err, "illegal character token '%s' on line %lu in file "
                     "'%s'", gt_str_get(token),
                     gt_tokenizer_get_line_number(tz),
                     gt_tokenizer_get_filename(tz));
        had_err = -1;
        break;
      }
      gt_str_delete(token);
      gt_tokenizer_next_token(tz);
      gt_assert(!had_err);
      return 0;
    }
    gt_str_delete(token);
    gt_tokenizer_next_token(tz);
  }
  if (!had_err) {
    if (!gt_array_size(index_to_alpha_char_mapping)) {
      gt_error_set(err, "could not parse a single alphabet character in file "
                   "'%s' (file empty or directory?)",
                   gt_tokenizer_get_filename(tz));
    had_err = -1;
    }
  }
  gt_str_delete(token);
  return had_err;
}

static int parse_score_line(GtScoreMatrix *sm, GtTokenizer *tz,
                            GtArray *index_to_alpha_char_mapping,
                            char *parsed_characters, GtError *err)
{
  unsigned int num_of_chars, i = 0;
  char amino_acid;
  int score, had_err = 0;
  GtStr *token;
  gt_assert(sm && tz && index_to_alpha_char_mapping);
  gt_error_check(err);
  token = gt_tokenizer_get_token(tz);
  gt_assert(token);
  if (gt_str_length(token) != 1) {
    gt_error_set(err, "illegal character token '%s' on line %lu in file '%s'",
                 gt_str_get(token), gt_tokenizer_get_line_number(tz),
                 gt_tokenizer_get_filename(tz));
    had_err = -1;
  }
  amino_acid = gt_str_get(token)[0];
  /* check for character duplications */
  if (parsed_characters[(int) amino_acid]) {
    gt_error_set(err, "multiple character '%c' entry on line %lu in file '%s'",
                 amino_acid, gt_tokenizer_get_line_number(tz),
                 gt_tokenizer_get_filename(tz));
    had_err = -1;
  }
  parsed_characters[(int) amino_acid] = GT_UNDEF_CHAR;
  gt_str_delete(token);
  if (!had_err) {
    num_of_chars = gt_alphabet_num_of_chars(sm->alphabet);
    gt_tokenizer_next_token(tz);
    while ((token = gt_tokenizer_get_token(tz))) {
      unsigned int idx1, idx2;
      /* the tokenizer can return tokens which are empty except for a newline
         -> skip these */
      if (!strcmp(gt_str_get(token), "\n")) {
        gt_str_delete(token);
        gt_tokenizer_next_token(tz);
        if (gt_tokenizer_line_start(tz))
          break;
        continue;
      }
      /* token is not empty -> parse score */
      had_err = gt_parse_int_line(&score, gt_str_get(token),
                                  gt_tokenizer_get_line_number(tz),
                                  gt_tokenizer_get_filename(tz), err);
      if (had_err)
        break;
      idx1 = gt_alphabet_encode(sm->alphabet, amino_acid);
      idx2 = gt_alphabet_encode(sm->alphabet, *(char*)
                                gt_array_get(index_to_alpha_char_mapping, i));
      gt_score_matrix_set_score(sm,
                                idx1 == WILDCARD ? num_of_chars : idx1,
                                idx2 == WILDCARD ? num_of_chars : idx2,
                                score);
      i++;
      gt_str_delete(token);
      gt_tokenizer_next_token(tz);
      if (gt_tokenizer_line_start(tz))
        break;
    }
  }
  return had_err;
}

/* the score matrix parser */
static int parse_score_matrix(GtScoreMatrix *sm, const char *path,
                              GtError *err)
{
  GtTokenizer *tz;
  GtArray *index_to_alpha_char_mapping;
  unsigned int parsed_score_lines = 0;
  char parsed_characters[UCHAR_MAX] = { 0 };
  int had_err = 0;
  gt_error_check(err);
  gt_assert(sm && path && sm->alphabet);
  tz = gt_tokenizer_new(gt_io_new(path, "r"));
  index_to_alpha_char_mapping = gt_array_new(sizeof (char));
  gt_tokenizer_skip_comment_lines(tz);
  had_err = parse_alphabet_line(index_to_alpha_char_mapping, tz, err);
  if (!had_err) {
    while (gt_tokenizer_has_token(tz)) {
      had_err = parse_score_line(sm, tz, index_to_alpha_char_mapping,
                                 parsed_characters, err);
      if (had_err)
        break;
      parsed_score_lines++;
    }
  }

  /* check the number of parsed score lines */
  if (!had_err &&
      parsed_score_lines != gt_array_size(index_to_alpha_char_mapping)) {
    gt_error_set(err, "the score matrix given in '%s' is not symmetric", path);
    had_err = -1;
  }

  gt_array_delete(index_to_alpha_char_mapping);
  gt_tokenizer_delete(tz);

  return had_err;
}

GtScoreMatrix* gt_score_matrix_new_read_protein(const char *path, GtError *err)
{
  GtAlphabet *protein_alpha;
  GtScoreMatrix *sm;
  int had_err;

  gt_error_check(err);
  gt_assert(path);

  /* create score matrix */
  protein_alpha = gt_alphabet_new_protein();
  sm = gt_score_matrix_new(protein_alpha);
  gt_alphabet_delete(protein_alpha);

  /* parse matrix file */
  had_err = parse_score_matrix(sm, path, err);

  if (had_err) {
    gt_score_matrix_delete(sm);
    return NULL;
  }
  return sm;
}

GtScoreMatrix* gt_score_matrix_new_read(const char *path, GtAlphabet *alphabet,
                                        GtError *err)
{
  GtScoreMatrix *sm;
  gt_error_check(err);
  gt_assert(path && alphabet);
  sm = gt_score_matrix_new(alphabet);
  if (parse_score_matrix(sm, path, err)) {
    gt_score_matrix_delete(sm);
    return NULL;
  }
  return sm;
}

unsigned int gt_score_matrix_get_dimension(const GtScoreMatrix *sm)
{
  gt_assert(sm);
  return sm->dimension;
}

int gt_score_matrix_get_score(const GtScoreMatrix *sm,
                              unsigned int idx1, unsigned int idx2)
{
  gt_assert(sm);
  idx1 = (idx1 == WILDCARD) ? sm->dimension - 1 : idx1;
  idx2 = (idx2 == WILDCARD) ? sm->dimension - 1 : idx2;
  /* indices are valid */
  gt_assert(idx1 < sm->dimension && idx2 < sm->dimension);
  return sm->scores[idx1][idx2];
}

void gt_score_matrix_set_score(GtScoreMatrix *sm,
                               unsigned int idx1, unsigned int idx2, int score)
{
  gt_assert(sm);
  idx1 = (idx1 == WILDCARD) ? sm->dimension - 1 : idx1;
  idx2 = (idx2 == WILDCARD) ? sm->dimension - 1 : idx2;
  /* indices are valid */
  gt_assert(idx1 < sm->dimension && idx2 < sm->dimension);
  sm->scores[idx1][idx2] = score;
}

const int** gt_score_matrix_get_scores(const GtScoreMatrix *sm)
{
  gt_assert(sm);
  return (const int**) sm->scores;
}

void gt_score_matrix_show(const GtScoreMatrix *sm, FILE *fp)
{
  unsigned i, j;
  gt_assert(sm && fp);
  /* show alphabet line */
  gt_xfputc(' ', fp);
  for (i = 0; i < gt_alphabet_size(sm->alphabet); i++)
    fprintf(fp, "  %c", gt_alphabet_decode(sm->alphabet, i));
  gt_xfputc('\n', fp);
  /* show score lines */
  for (i = 0; i < gt_alphabet_size(sm->alphabet); i++) {
    gt_xfputc(gt_alphabet_decode(sm->alphabet, i), fp);
    for (j = 0; j < gt_alphabet_size(sm->alphabet); j++)
      fprintf(fp, " %2d", gt_score_matrix_get_score(sm, i, j));
    gt_xfputc('\n', fp);
  }
}

void gt_score_matrix_delete(GtScoreMatrix *sm)
{
  if (!sm) return;
  gt_alphabet_delete(sm->alphabet);
  gt_array2dim_delete(sm->scores);
  gt_free(sm);
}
