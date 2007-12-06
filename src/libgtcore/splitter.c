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
#include <stdio.h>
#include <string.h>
#include "libgtcore/dynalloc.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/splitter.h"
#include "libgtcore/xansi.h"

struct Splitter {
  char **tokens;
  unsigned long num_of_tokens;
  size_t allocated;
};

Splitter* splitter_new(void)
{
  return ma_calloc(1, sizeof (Splitter));
}

void splitter_split(Splitter *s, char *string, unsigned long length,
                    char delimiter)
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
                           (s->num_of_tokens + 1) * sizeof (char*));
    s->tokens[s->num_of_tokens++] = string_index;
    string_index = end_of_token + 1;
  }

  /* save last token */
  if ((s->num_of_tokens + 2) * sizeof (char*) > s->allocated)
    s->tokens = dynalloc(s->tokens, &s->allocated,
                         (s->num_of_tokens + 2) * sizeof (char*));
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

int splitter_unit_test(Error *err)
{
  static char string_1[]  = "a bb ccc dddd eeeee",
              string_2[]  = "a\tbb\tccc\tdddd\teeeee",
              string_3[]  = "",
              string_4[]  = "a  b",
              string_5[]  = "ac bc ",
              string_6[]  = "test";
  Splitter *s;
  int had_err = 0;
  error_check(err);
  s = splitter_new();

  /* string_1 */
  ensure(had_err, !splitter_size(s));
  splitter_split(s, string_1, strlen(string_1), ' ');
  ensure(had_err, splitter_size(s) == 5);
  ensure(had_err, strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 1), "bb") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 2), "ccc") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 3), "dddd") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 4), "eeeee") == 0);
  splitter_reset(s);

  /* string_2 */
  ensure(had_err, !splitter_size(s));
  splitter_split(s, string_2, strlen(string_2), '\t');
  ensure(had_err, splitter_size(s) == 5);
  ensure(had_err, strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 1), "bb") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 2), "ccc") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 3), "dddd") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 4), "eeeee") == 0);
  splitter_reset(s);

  /* string_3 */
  ensure(had_err, !splitter_size(s));
  splitter_split(s, string_3, strlen(string_3), '\t');
  ensure(had_err, splitter_size(s) == 1);
  ensure(had_err, strcmp(splitter_get_token(s, 0), "") == 0);
  splitter_reset(s);

  /* string_4 */
  ensure(had_err, !splitter_size(s));
  splitter_split(s, string_4, strlen(string_4), ' ');
  ensure(had_err, splitter_size(s) == 3);
  ensure(had_err, strcmp(splitter_get_token(s, 0), "a") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 1), "") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 2), "b") == 0);
  splitter_reset(s);

  /* string_5 */
  ensure(had_err, !splitter_size(s));
  splitter_split(s, string_5, strlen(string_5), ' ');
  ensure(had_err, splitter_size(s) == 3);
  ensure(had_err, strcmp(splitter_get_token(s, 0), "ac") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 1), "bc") == 0);
  ensure(had_err, strcmp(splitter_get_token(s, 2), "") == 0);
  splitter_reset(s);

  /* string_6 */
  ensure(had_err, !splitter_size(s));
  splitter_split(s, string_6, strlen(string_6), ';');
  ensure(had_err, splitter_size(s) == 1);
  ensure(had_err, strcmp(splitter_get_token(s, 0), "test") == 0);

  /* free */
  splitter_delete(s);

  return had_err;
}

void splitter_delete(Splitter *s)
{
  if (!s) return;
  ma_free(s->tokens);
  ma_free(s);
}
