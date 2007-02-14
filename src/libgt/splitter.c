/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "dynalloc.h"
#include "ensure.h"
#include "splitter.h"
#include "xansi.h"

struct Splitter {
  char **tokens;
  unsigned long num_of_tokens;
  size_t allocated;
};

Splitter* splitter_new(void)
{
  return xcalloc(1, sizeof(Splitter));
}

void splitter_split(Splitter *t, char *string, unsigned long length,
                    char delimiter)
{

  char *end_of_token, *string_index = string;

  assert(t && string);

  /* splitting */
  while (string_index < string + length &&
         (end_of_token = strchr(string_index, delimiter))) {
    assert(end_of_token != NULL);
    *end_of_token = '\0';
    if ((t->num_of_tokens + 1) * sizeof(char*) > t->allocated)
      t->tokens = dynalloc(t->tokens, &t->allocated,
                           (t->num_of_tokens + 1) * sizeof(char*));
    t->tokens[t->num_of_tokens++] = string_index;
    string_index = end_of_token + 1;
  }

  /* save last token */
  if ((t->num_of_tokens + 2) * sizeof(char*) > t->allocated)
    t->tokens = dynalloc(t->tokens, &t->allocated,
                         (t->num_of_tokens + 2) * sizeof(char*));
  t->tokens[t->num_of_tokens++] = string_index;
  t->tokens[t->num_of_tokens]   = NULL;

  assert(t->num_of_tokens);
}

char** splitter_get_tokens(Splitter *t)
{
  assert(t);
  return t->tokens;
}

char* splitter_get_token(Splitter *t, unsigned long token_num)
{
  assert(t && token_num < t->num_of_tokens);
  return t->tokens[token_num];
}

void splitter_reset(Splitter *t)
{
  assert(t);
  if (t->tokens) t->tokens[0] = NULL;
  t->num_of_tokens = 0;
}

unsigned long splitter_size(Splitter *t)
{
  assert(t);
  return t->num_of_tokens;
}

int splitter_unit_test(void)
{
  Splitter *s = splitter_new();
  static char string_1[]  = "a bb ccc dddd eeeee",
              string_2[]  = "a\tbb\tccc\tdddd\teeeee",
              string_3[]  = "",
              string_4[]  = "a  b",
              string_5[]  = "ac bc ",
              string_6[]  = "test";

  /* string_1 */
  ensure(!splitter_size(s));
  splitter_split(s, string_1, strlen(string_1), ' ');
  ensure(splitter_size(s) == 5);
  ensure(strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(strcmp(splitter_get_token(s, 1), "bb") == 0);
  ensure(strcmp(splitter_get_token(s, 2), "ccc") == 0);
  ensure(strcmp(splitter_get_token(s, 3), "dddd") == 0);
  ensure(strcmp(splitter_get_token(s, 4), "eeeee") == 0);
  splitter_reset(s);

  /* string_2 */
  ensure(!splitter_size(s));
  splitter_split(s, string_2, strlen(string_2), '\t');
  ensure(splitter_size(s) == 5);
  ensure(strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(strcmp(splitter_get_token(s, 1), "bb") == 0);
  ensure(strcmp(splitter_get_token(s, 2), "ccc") == 0);
  ensure(strcmp(splitter_get_token(s, 3), "dddd") == 0);
  ensure(strcmp(splitter_get_token(s, 4), "eeeee") == 0);
  splitter_reset(s);

  /* string_3 */
  ensure(!splitter_size(s));
  splitter_split(s, string_3, strlen(string_3), '\t');
  ensure(splitter_size(s) == 1);
  ensure(strcmp(splitter_get_token(s, 0), "") == 0);
  splitter_reset(s);

  /* string_4 */
  ensure(!splitter_size(s));
  splitter_split(s, string_4, strlen(string_4), ' ');
  ensure(splitter_size(s) == 3);
  ensure(strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(strcmp(splitter_get_token(s, 1), "") == 0);
  ensure(strcmp(splitter_get_token(s, 2), "b") == 0);
  splitter_reset(s);

  /* string_5 */
  ensure(!splitter_size(s));
  splitter_split(s, string_5, strlen(string_5), ' ');
  ensure(splitter_size(s) == 3);
  ensure(strcmp(splitter_get_token(s, 0), "ac") == 0);
  ensure(strcmp(splitter_get_token(s, 1), "bc") == 0);
  ensure(strcmp(splitter_get_token(s, 2), "") == 0);
  splitter_reset(s);

  /* string_6 */
  ensure(!splitter_size(s));
  splitter_split(s, string_6, strlen(string_6), ';');
  ensure(splitter_size(s) == 1);
  ensure(strcmp(splitter_get_token(s, 0), "test") == 0);

  /* free */
  splitter_free(s);

  return EXIT_SUCCESS;
}

void splitter_free(Splitter *t)
{
  if (!t) return;
  free(t->tokens);
  free(t);
}
