/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "array.h"
#include "array2dim.h"
#include "parseutils.h"
#include "scorematrix.h"
#include "str.h"
#include "tokenizer.h"
#include "undef.h"
#include "xansi.h"

struct ScoreMatrix {
  Alpha *alpha;
  unsigned int dimension;
  int **scores;
};

ScoreMatrix* scorematrix_new(Alpha *alpha)
{
  ScoreMatrix *s;
  assert(alpha);
  s = xmalloc(sizeof (ScoreMatrix));
  s->alpha = alpha_ref(alpha);
  s->dimension = alpha_size(alpha);
  array2dim_calloc(s->scores, s->dimension, s->dimension, int);
  return s;
}

static int parse_alphabet_line(Array *index_to_alpha_char_mapping,
                                Tokenizer *tz, Error *err)
{
  Str *token;
  char *tokenstr, amino_acid, parsed_characters[UCHAR_MAX] = { 0 };
  int has_err = 0;
  error_check(err);
  assert(index_to_alpha_char_mapping && tz);
  assert(!array_size(index_to_alpha_char_mapping));
  while ((token = tokenizer_get_token(tz))) {
    if (str_length(token) > 2) {
      error_set(err, "illegal character token '%s' on line %lu in file '%s'",
                str_get(token), tokenizer_get_line_number(tz),
                tokenizer_get_filename(tz));
      has_err = -1;
      break;
    }
    tokenstr = str_get(token);
    amino_acid = tokenstr[0];
    /* check for character duplications */
    if (parsed_characters[(int) amino_acid]) {
      error_set(err, "the character '%c' appears more then once on line %lu in "
                "file  '%s'", amino_acid, tokenizer_get_line_number(tz),
                tokenizer_get_filename(tz));
      has_err = -1;
      break;
    }
    parsed_characters[(int) amino_acid] = UNDEFCHAR;
    if (amino_acid == '\n') {
      str_delete(token);
      tokenizer_next_token(tz);
      assert(!has_err);
      return 0;
    }
    array_add(index_to_alpha_char_mapping, amino_acid);
    if (str_length(token) == 2) {
      if (tokenstr[1] != '\n') {
        error_set(err, "illegal character token '%s' on line %lu in file '%s'",
                  str_get(token), tokenizer_get_line_number(tz),
                  tokenizer_get_filename(tz));
        has_err = -1;
        break;
      }
      str_delete(token);
      tokenizer_next_token(tz);
      assert(!has_err);
      return 0;
    }
    str_delete(token);
    tokenizer_next_token(tz);
  }
  if (!has_err) {
    if (!array_size(index_to_alpha_char_mapping)) {
      error_set(err, "could not parse a single alphabet character in file '%s' "
                "(file empty or directory?)", tokenizer_get_filename(tz));
    has_err = -1;
    }
  }
  str_delete(token);
  return has_err;
}

static int parse_score_line(ScoreMatrix *s, Tokenizer *tz,
                            Array *index_to_alpha_char_mapping,
                            char *parsed_characters, Error *err)
{
  unsigned int i = 0;
  char amino_acid;
  int score, has_err = 0;
  Str *token;
  assert(s && tz && index_to_alpha_char_mapping);
  error_check(err);
  token = tokenizer_get_token(tz);
  assert(token);
  if (str_length(token) != 1) {
    error_set(err, "illegal character token '%s' on line %lu in file '%s'",
              str_get(token), tokenizer_get_line_number(tz),
              tokenizer_get_filename(tz));
    has_err = -1;
  }
  amino_acid = str_get(token)[0];
  /* check for character duplications */
  if (parsed_characters[(int) amino_acid]) {
    error_set(err, "multiple character '%c' entry on line %lu in file '%s'",
              amino_acid, tokenizer_get_line_number(tz),
              tokenizer_get_filename(tz));
    has_err = -1;
  }
  parsed_characters[(int) amino_acid] = UNDEFCHAR;
  str_delete(token);
  if (!has_err) {
    tokenizer_next_token(tz);
    while ((token = tokenizer_get_token(tz))) {
      has_err = parse_int(&score, str_get(token), tokenizer_get_line_number(tz),
                          tokenizer_get_filename(tz), err);
      if (has_err)
        break;
      scorematrix_set_score(s,
                            (unsigned char) alpha_encode(s->alpha, amino_acid),
                            (unsigned char) alpha_encode(s->alpha, *(char*)
                            array_get(index_to_alpha_char_mapping, i)), score);
      i++;
      str_delete(token);
      tokenizer_next_token(tz);
      if (tokenizer_line_start(tz))
          break;
    }
  }
  return has_err;
}

/* the score matrix parser */
static int parse_scorematrix(ScoreMatrix *s, const char *path, Error *err)
{
  Tokenizer *tz;
  Array *index_to_alpha_char_mapping;
  unsigned int parsed_score_lines = 0;
  char parsed_characters[UCHAR_MAX] = { 0 };
  int has_err = 0;
  assert(s && path && s->alpha);
  error_check(err);
  tz = tokenizer_new(io_new(path, "r"));
  index_to_alpha_char_mapping = array_new(sizeof (char));
  tokenizer_skip_comment_lines(tz);
  has_err = parse_alphabet_line(index_to_alpha_char_mapping, tz, err);
  if (!has_err) {
    while (tokenizer_has_token(tz)) {
      has_err = parse_score_line(s, tz, index_to_alpha_char_mapping,
                                 parsed_characters, err);
      if (has_err)
        break;
      parsed_score_lines++;
    }
  }

  /* check the number of parsed score lines */
  if (!has_err &&
      parsed_score_lines != array_size(index_to_alpha_char_mapping)) {
    error_set(err, "the scorematrix given in '%s' is not symmetric", path);
    has_err = -1;
  }

  array_delete(index_to_alpha_char_mapping);
  tokenizer_delete(tz);

  return has_err;
}

ScoreMatrix* scorematrix_read_protein(const char *path, Error *err)
{
  Alpha *protein_alpha;
  ScoreMatrix *s;
  int has_err;
  assert(path);
  error_check(err);

  /* create score matrix */
  protein_alpha = alpha_new_protein();
  s = scorematrix_new(protein_alpha);
  alpha_delete(protein_alpha);

  /* parse matrix file */
  has_err = parse_scorematrix(s, path, err);

  if (has_err) {
    scorematrix_delete(s);
    return NULL;
  }
  return s;
}

int scorematrix_get_score(const ScoreMatrix *s,
                          unsigned char idx1, unsigned char idx2)
{
  assert(s);
  assert(idx1 < s->dimension && idx2 < s->dimension); /* indices are valid */
  return s->scores[idx1][idx2];
}

void scorematrix_set_score(ScoreMatrix *s,
                           unsigned char idx1, unsigned char idx2, int score)
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
  array2dim_delete(s->scores);
  free(s);
}
