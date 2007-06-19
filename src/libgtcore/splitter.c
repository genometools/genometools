/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <libgtcore/dynalloc.h>
#include <libgtcore/ensure.h>
#include <libgtcore/splitter.h>
#include <libgtcore/xansi.h>

struct Splitter {
  char **tokens;
  unsigned long num_of_tokens;
  size_t allocated;
};

Splitter* splitter_new(Env *env)
{
  return env_ma_calloc(env, 1, sizeof (Splitter));
}

void splitter_split(Splitter *s, char *string, unsigned long length,
                    char delimiter, Env *env)
{

  char *end_of_token, *string_index = string;

  assert(s && string);

  /* splitting */
  while (string_index < string + length &&
         (end_of_token = strchr(string_index, delimiter))) {
    assert(end_of_token);
    *end_of_token = '\0';
    if ((s->num_of_tokens + 1) * sizeof (char*) > s->allocated)
      s->tokens = dynalloc(s->tokens, &s->allocated,
                           (s->num_of_tokens + 1) * sizeof (char*), env);
    s->tokens[s->num_of_tokens++] = string_index;
    string_index = end_of_token + 1;
  }

  /* save last token */
  if ((s->num_of_tokens + 2) * sizeof (char*) > s->allocated)
    s->tokens = dynalloc(s->tokens, &s->allocated,
                         (s->num_of_tokens + 2) * sizeof (char*), env);
  s->tokens[s->num_of_tokens++] = string_index;
  s->tokens[s->num_of_tokens]   = NULL;

  assert(s->num_of_tokens);
}

char** splitter_get_tokens(Splitter *s)
{
  assert(s);
  return s->tokens;
}

char* splitter_get_token(Splitter *s, unsigned long token_num)
{
  assert(s && token_num < s->num_of_tokens);
  return s->tokens[token_num];
}

void splitter_reset(Splitter *s)
{
  assert(s);
  if (s->tokens) s->tokens[0] = NULL;
  s->num_of_tokens = 0;
}

unsigned long splitter_size(Splitter *s)
{
  assert(s);
  return s->num_of_tokens;
}

int splitter_unit_test(Env *env)
{
  static char string_1[]  = "a bb ccc dddd eeeee",
              string_2[]  = "a\tbb\tccc\tdddd\teeeee",
              string_3[]  = "",
              string_4[]  = "a  b",
              string_5[]  = "ac bc ",
              string_6[]  = "test";
  Splitter *s;
  int has_err = 0;
  env_error_check(env);
  s = splitter_new(env);

  /* string_1 */
  ensure(has_err, !splitter_size(s));
  splitter_split(s, string_1, strlen(string_1), ' ', env);
  ensure(has_err, splitter_size(s) == 5);
  ensure(has_err, strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 1), "bb") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 2), "ccc") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 3), "dddd") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 4), "eeeee") == 0);
  splitter_reset(s);

  /* string_2 */
  ensure(has_err, !splitter_size(s));
  splitter_split(s, string_2, strlen(string_2), '\t', env);
  ensure(has_err, splitter_size(s) == 5);
  ensure(has_err, strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 1), "bb") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 2), "ccc") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 3), "dddd") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 4), "eeeee") == 0);
  splitter_reset(s);

  /* string_3 */
  ensure(has_err, !splitter_size(s));
  splitter_split(s, string_3, strlen(string_3), '\t', env);
  ensure(has_err, splitter_size(s) == 1);
  ensure(has_err, strcmp(splitter_get_token(s, 0), "") == 0);
  splitter_reset(s);

  /* string_4 */
  ensure(has_err, !splitter_size(s));
  splitter_split(s, string_4, strlen(string_4), ' ', env);
  ensure(has_err, splitter_size(s) == 3);
  ensure(has_err, strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 1), "") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 2), "b") == 0);
  splitter_reset(s);

  /* string_5 */
  ensure(has_err, !splitter_size(s));
  splitter_split(s, string_5, strlen(string_5), ' ', env);
  ensure(has_err, splitter_size(s) == 3);
  ensure(has_err, strcmp(splitter_get_token(s, 0), "ac") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 1), "bc") == 0);
  ensure(has_err, strcmp(splitter_get_token(s, 2), "") == 0);
  splitter_reset(s);

  /* string_6 */
  ensure(has_err, !splitter_size(s));
  splitter_split(s, string_6, strlen(string_6), ';', env);
  ensure(has_err, splitter_size(s) == 1);
  ensure(has_err, strcmp(splitter_get_token(s, 0), "test") == 0);

  /* free */
  splitter_delete(s, env);

  return has_err;
}

void splitter_delete(Splitter *s, Env *env)
{
  if (!s) return;
  env_ma_free(s->tokens, env);
  env_ma_free(s, env);
}
