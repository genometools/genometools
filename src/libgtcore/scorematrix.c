/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/array.h"
#include "libgtcore/array2dim.h"
#include "libgtcore/parseutils.h"
#include "libgtcore/scorematrix.h"
#include "libgtcore/str.h"
#include "libgtcore/tokenizer.h"
#include "libgtcore/undef.h"
#include "libgtcore/xansi.h"

struct ScoreMatrix {
  Alpha *alpha;
  unsigned int dimension;
  int **scores;
};

ScoreMatrix* scorematrix_new(Alpha *alpha)
{
  ScoreMatrix *s;
  assert(alpha);
  s = ma_malloc(sizeof (ScoreMatrix));
  s->alpha = alpha_ref(alpha);
  s->dimension = alpha_size(alpha);
  array2dim_calloc(s->scores, s->dimension, s->dimension);
  return s;
}

static int parse_alphabet_line(Array *index_to_alpha_char_mapping,
                               Tokenizer *tz, Error *e)
{
  Str *token;
  char *tokenstr, amino_acid, parsed_characters[UCHAR_MAX] = { 0 };
  int had_err = 0;
  error_check(e);
  assert(index_to_alpha_char_mapping && tz);
  assert(!array_size(index_to_alpha_char_mapping));
  while ((token = tokenizer_get_token(tz))) {
    if (str_length(token) > 2) {
      error_set(e, "illegal character token '%s' on line %lu in file '%s'",
                str_get(token), tokenizer_get_line_number(tz),
                tokenizer_get_filename(tz));
      had_err = -1;
      break;
    }
    tokenstr = str_get(token);
    amino_acid = tokenstr[0];
    /* check for character duplications */
    if (parsed_characters[(int) amino_acid]) {
      error_set(e, "the character '%c' appears more then once on line %lu in "
                "file  '%s'", amino_acid, tokenizer_get_line_number(tz),
                tokenizer_get_filename(tz));
      had_err = -1;
      break;
    }
    parsed_characters[(int) amino_acid] = UNDEF_CHAR;
    if (amino_acid == '\n') {
      str_delete(token);
      tokenizer_next_token(tz);
      assert(!had_err);
      return 0;
    }
    array_add(index_to_alpha_char_mapping, amino_acid);
    if (str_length(token) == 2) {
      if (tokenstr[1] != '\n') {
        error_set(e, "illegal character token '%s' on line %lu in file '%s'",
                  str_get(token), tokenizer_get_line_number(tz),
                  tokenizer_get_filename(tz));
        had_err = -1;
        break;
      }
      str_delete(token);
      tokenizer_next_token(tz);
      assert(!had_err);
      return 0;
    }
    str_delete(token);
    tokenizer_next_token(tz);
  }
  if (!had_err) {
    if (!array_size(index_to_alpha_char_mapping)) {
      error_set(e, "could not parse a single alphabet character in file "
                "'%s' (file empty or directory?)", tokenizer_get_filename(tz));
    had_err = -1;
    }
  }
  str_delete(token);
  return had_err;
}

static int parse_score_line(ScoreMatrix *s, Tokenizer *tz,
                            Array *index_to_alpha_char_mapping,
                            char *parsed_characters, Error *e)
{
  unsigned int i = 0;
  char amino_acid;
  int score, had_err = 0;
  Str *token;
  assert(s && tz && index_to_alpha_char_mapping);
  error_check(e);
  token = tokenizer_get_token(tz);
  assert(token);
  if (str_length(token) != 1) {
    error_set(e, "illegal character token '%s' on line %lu in file '%s'",
              str_get(token), tokenizer_get_line_number(tz),
              tokenizer_get_filename(tz));
    had_err = -1;
  }
  amino_acid = str_get(token)[0];
  /* check for character duplications */
  if (parsed_characters[(int) amino_acid]) {
    error_set(e, "multiple character '%c' entry on line %lu in file '%s'",
              amino_acid, tokenizer_get_line_number(tz),
              tokenizer_get_filename(tz));
    had_err = -1;
  }
  parsed_characters[(int) amino_acid] = UNDEF_CHAR;
  str_delete(token);
  if (!had_err) {
    tokenizer_next_token(tz);
    while ((token = tokenizer_get_token(tz))) {
      had_err = parse_int_line(&score, str_get(token),
                               tokenizer_get_line_number(tz),
                               tokenizer_get_filename(tz), e);
      if (had_err)
        break;
      scorematrix_set_score(s,
                            alpha_encode(s->alpha, amino_acid),
                            alpha_encode(s->alpha, *(char*)
                            array_get(index_to_alpha_char_mapping, i)), score);
      i++;
      str_delete(token);
      tokenizer_next_token(tz);
      if (tokenizer_line_start(tz))
          break;
    }
  }
  return had_err;
}

/* the score matrix parser */
static int parse_scorematrix(ScoreMatrix *s, const char *path, Error *e)
{
  Tokenizer *tz;
  Array *index_to_alpha_char_mapping;
  unsigned int parsed_score_lines = 0;
  char parsed_characters[UCHAR_MAX] = { 0 };
  int had_err = 0;
  error_check(e);
  assert(s && path && s->alpha);
  tz = tokenizer_new(io_new(path, "r"));
  index_to_alpha_char_mapping = array_new(sizeof (char));
  tokenizer_skip_comment_lines(tz);
  had_err = parse_alphabet_line(index_to_alpha_char_mapping, tz, e);
  if (!had_err) {
    while (tokenizer_has_token(tz)) {
      had_err = parse_score_line(s, tz, index_to_alpha_char_mapping,
                                 parsed_characters, e);
      if (had_err)
        break;
      parsed_score_lines++;
    }
  }

  /* check the number of parsed score lines */
  if (!had_err &&
      parsed_score_lines != array_size(index_to_alpha_char_mapping)) {
    error_set(e, "the scorematrix given in '%s' is not symmetric", path);
    had_err = -1;
  }

  array_delete(index_to_alpha_char_mapping);
  tokenizer_delete(tz);

  return had_err;
}

ScoreMatrix* scorematrix_read_protein(const char *path, Error *e)
{
  Alpha *protein_alpha;
  ScoreMatrix *s;
  int had_err;

  error_check(e);
  assert(path);

  /* create score matrix */
  protein_alpha = alpha_new_protein();
  s = scorematrix_new(protein_alpha);
  alpha_delete(protein_alpha);

  /* parse matrix file */
  had_err = parse_scorematrix(s, path, e);

  if (had_err) {
    scorematrix_delete(s);
    return NULL;
  }
  return s;
}

unsigned int scorematrix_get_dimension(const ScoreMatrix *s)
{
  assert(s);
  return s->dimension;
}

int scorematrix_get_score(const ScoreMatrix *s,
                          unsigned int idx1, unsigned int idx2)
{
  assert(s);
  assert(idx1 < s->dimension && idx2 < s->dimension); /* indices are valid */
  return s->scores[idx1][idx2];
}

void scorematrix_set_score(ScoreMatrix *s,
                           unsigned int idx1, unsigned int idx2, int score)
{
  assert(s);
  assert(idx1 < s->dimension && idx2 < s->dimension); /* indices are valid */
  s->scores[idx1][idx2] = score;
}

void scorematrix_show(const ScoreMatrix *s, FILE *fp)
{
  unsigned i, j;
  assert(s && fp);
  /* show alphabet line */
  xfputc(' ', fp);
  for (i = 0; i < alpha_size(s->alpha); i++)
    fprintf(fp, "  %c", alpha_decode(s->alpha, i));
  xfputc('\n', fp);
  /* show score lines */
  for (i = 0; i < alpha_size(s->alpha); i++) {
    xfputc(alpha_decode(s->alpha, i), fp);
    for (j = 0; j < alpha_size(s->alpha); j++)
      fprintf(fp, " %2d", scorematrix_get_score(s, i, j));
    xfputc('\n', fp);
  }
}

void scorematrix_delete(ScoreMatrix *s)
{
  if (!s) return;
  alpha_delete(s->alpha);
  array2dim_delete(s->scores);
  ma_free(s);
}
