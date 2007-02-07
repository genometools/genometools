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
  s = xmalloc(sizeof(ScoreMatrix));
  s->alpha = alpha_ref(alpha);
  s->dimension = alpha_size(alpha);
  array2dim_calloc(s->scores, s->dimension, s->dimension, int);
  return s;
}

static void parse_alphabet_line(Array *index_to_alpha_char_mapping,
                                Tokenizer *tz)
{
  Str *token;
  char *tokenstr, amino_acid, parsed_characters[UCHAR_MAX] = { 0 };
  assert(index_to_alpha_char_mapping && tz);
  assert(!array_size(index_to_alpha_char_mapping));
  while ((token = tokenizer_get_token(tz))) {
    if (str_length(token) > 2) {
      error("illegal character token '%s' on line %lu in file '%s'",
            str_get(token), tokenizer_get_line_number(tz),
            tokenizer_get_filename(tz));
    }
    tokenstr = str_get(token);
    amino_acid = tokenstr[0];
    /* check for character duplications */
    if (parsed_characters[(int) amino_acid]) {
      error("the character '%c' appears more then once on line %lu in file "
            " '%s'", amino_acid, tokenizer_get_line_number(tz),
            tokenizer_get_filename(tz));
   }
    parsed_characters[(int) amino_acid] = UNDEFCHAR;
    if (amino_acid == '\n') {
      str_free(token);
      tokenizer_next_token(tz);
      return;
    }
    array_add(index_to_alpha_char_mapping, amino_acid);
    if (str_length(token) == 2) {
      if (tokenstr[1] != '\n') {
        error("illegal character token '%s' on line %lu in file '%s'",
              str_get(token), tokenizer_get_line_number(tz),
              tokenizer_get_filename(tz));
      }
      str_free(token);
      tokenizer_next_token(tz);
      return;
    }
    str_free(token);
    tokenizer_next_token(tz);
  }
  if (!array_size(index_to_alpha_char_mapping)) {
    error("could not parse a single alphabet character in file '%s' (file "
          "empty or directory?)", tokenizer_get_filename(tz));
  }
}

static void parse_score_line(ScoreMatrix *s, Tokenizer *tz,
                             Array *index_to_alpha_char_mapping,
                             char *parsed_characters)
{
  unsigned int i = 0;
  char amino_acid;
  int score;
  Str *token;
  assert(s && tz && index_to_alpha_char_mapping);
  token = tokenizer_get_token(tz);
  assert(token);
  if (str_length(token) != 1) {
    error("illegal character token '%s' on line %lu in file '%s'",
          str_get(token), tokenizer_get_line_number(tz),
          tokenizer_get_filename(tz));
  }
  amino_acid = str_get(token)[0];
  /* check for character duplications */
  if (parsed_characters[(int) amino_acid]) {
    error("multiple character '%c' entry on line %lu in file '%s'",
          amino_acid, tokenizer_get_line_number(tz),
          tokenizer_get_filename(tz));
  }
  parsed_characters[(int) amino_acid] = UNDEFCHAR;
  str_free(token);
  tokenizer_next_token(tz);
  while ((token = tokenizer_get_token(tz))) {
    score = parse_int(str_get(token), tokenizer_get_line_number(tz),
                      tokenizer_get_filename(tz), NULL);
    scorematrix_set_score(s, (unsigned char) alpha_encode(s->alpha, amino_acid),
                          (unsigned char) alpha_encode(s->alpha, *(char*)
                                          array_get(index_to_alpha_char_mapping,
                                                    i)), score);
    i++;
    str_free(token);
    tokenizer_next_token(tz);
    if (tokenizer_line_start(tz))
      break;
  }
}

/* the score matrix parser */
static void parse_scorematrix(ScoreMatrix *s, const char *path)
{
  Tokenizer *tz;
  Array *index_to_alpha_char_mapping;
  unsigned int parsed_score_lines = 0;
  char parsed_characters[UCHAR_MAX] = { 0 };
  assert(s && path && s->alpha);
  tz = tokenizer_new(io_new(path, "r"));
  index_to_alpha_char_mapping = array_new(sizeof(char));
  tokenizer_skip_comment_lines(tz);
  parse_alphabet_line(index_to_alpha_char_mapping, tz);
  while (tokenizer_has_token(tz)) {
    parse_score_line(s, tz, index_to_alpha_char_mapping, parsed_characters);
    parsed_score_lines++;
  }
  /* check the number of parsed score lines */
  if (parsed_score_lines != array_size(index_to_alpha_char_mapping))
    error("the scorematrix given in '%s' is not symmetric", path);
  array_free(index_to_alpha_char_mapping);
  tokenizer_free(tz);
}

ScoreMatrix* scorematrix_read_protein(const char *path)
{
  Alpha *protein_alpha;
  ScoreMatrix *s;
  assert(path);

  /* create score matrix */
  protein_alpha = alpha_new_protein();
  s = scorematrix_new(protein_alpha);
  alpha_free(protein_alpha);

  /* parse matrix file */
  parse_scorematrix(s, path);

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

void scorematrix_free(ScoreMatrix *s)
{
  if (!s) return;
  array2dim_free(s->scores);
  free(s);
}
