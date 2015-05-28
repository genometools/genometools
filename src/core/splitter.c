/*
  Copyright (c) 2006-2007 Gordon Gremme <gordon@gremme.org>
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

#include <stdio.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/dynalloc.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/splitter.h"
#include "core/xansi_api.h"

struct GtSplitter {
  char **tokens;
  GtUword num_of_tokens;
  size_t allocated;
};

GtSplitter* gt_splitter_new(void)
{
  return gt_calloc(1, sizeof (GtSplitter));
}

void gt_splitter_split(GtSplitter *s, char *string, GtUword length,
                       char delimiter)
{

  char *end_of_token, *string_index = string;

  gt_assert(s && string);

  /* splitting */
  while (string_index < string + length &&
         (end_of_token = strchr(string_index, delimiter))) {
    gt_assert(end_of_token);
    *end_of_token = '\0';
    if ((s->num_of_tokens + 1) * sizeof (char*) > s->allocated)
      s->tokens = gt_dynalloc(s->tokens, &s->allocated,
                              (s->num_of_tokens + 1) * sizeof (char*));
    s->tokens[s->num_of_tokens++] = string_index;
    string_index = end_of_token + 1;
  }

  /* save last token */
  if ((s->num_of_tokens + 2) * sizeof (char*) > s->allocated)
    s->tokens = gt_dynalloc(s->tokens, &s->allocated,
                            (s->num_of_tokens + 2) * sizeof (char*));
  s->tokens[s->num_of_tokens++] = string_index;
  s->tokens[s->num_of_tokens]   = NULL;

  gt_assert(s->num_of_tokens);
}

void gt_splitter_split_non_empty(GtSplitter *s, char *string, GtUword length,
                                 char delimiter)
{

  char *end_of_token, *string_index = string;

  gt_assert(s && string);

  /* splitting */
  while (string_index < string + length &&
         (end_of_token = strchr(string_index, delimiter))) {
    gt_assert(end_of_token);
    *end_of_token = '\0';
    /* skip empty token */
    if (end_of_token > string_index) {
      if ((s->num_of_tokens + 2) * sizeof (*s->tokens) > s->allocated)
        s->tokens = gt_dynalloc(s->tokens, &s->allocated,
                                (s->num_of_tokens + 2) * sizeof (*s->tokens));
      s->tokens[s->num_of_tokens++] = string_index;
    }
    string_index = end_of_token + 1;
  }

  /* save last token if not empty*/
  if (string_index < string + length) {
    if ((s->num_of_tokens + 2) * sizeof (*s->tokens) > s->allocated)
      s->tokens = gt_dynalloc(s->tokens, &s->allocated,
                              (s->num_of_tokens + 2) * sizeof (*s->tokens));
    s->tokens[s->num_of_tokens++] = string_index;
  }
  s->tokens[s->num_of_tokens] = NULL;
}

char** gt_splitter_get_tokens(GtSplitter *s)
{
  gt_assert(s);
  return s->tokens;
}

char* gt_splitter_get_token(GtSplitter *s, GtUword token_num)
{
  gt_assert(s && token_num < s->num_of_tokens);
  return s->tokens[token_num];
}

void gt_splitter_reset(GtSplitter *s)
{
  gt_assert(s);
  if (s->tokens) s->tokens[0] = NULL;
  s->num_of_tokens = 0;
}

GtUword gt_splitter_size(GtSplitter *s)
{
  gt_assert(s);
  return s->num_of_tokens;
}

void gt_splitter_delete(GtSplitter *s)
{
  if (!s) return;
  gt_free(s->tokens);
  gt_free(s);
}

#define RESET_DELIMITER(STR, DELIM) \
  STR[strlen(STR)] = DELIM

int gt_splitter_unit_test(GtError *err)
{
  static char string_1[]  = "a bb ccc dddd eeeee",
              string_2[]  = "a\tbb\tccc\tdddd\teeeee",
              string_3[]  = "",
              string_4[]  = "a  b",
              string_5[]  = "ac bc ",
              string_6[]  = "test",
              delimiter1 = ' ',
              delimiter2 = '\t',
              delimiter3 = ';';
  GtSplitter *s;
  int had_err = 0;
  GtUword idx;
  gt_error_check(err);
  s = gt_splitter_new();

  /* string_1 */
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split(s, string_1, strlen(string_1), delimiter1);
  gt_ensure(gt_splitter_size(s) == 5);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "a") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 1), "bb") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 2), "ccc") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 3), "dddd") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 4), "eeeee") == 0);
  for (idx = 0; idx < gt_splitter_size(s) - 1; idx++)
    RESET_DELIMITER(gt_splitter_get_token(s, idx), delimiter1);
  gt_splitter_reset(s);
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split_non_empty(s, string_1, strlen(string_1), delimiter1);
  gt_ensure(gt_splitter_size(s) == 5);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "a") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 1), "bb") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 2), "ccc") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 3), "dddd") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 4), "eeeee") == 0);
  gt_splitter_reset(s);

  /* string_2 */
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split(s, string_2, strlen(string_2), delimiter2);
  gt_ensure(gt_splitter_size(s) == 5);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "a") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 1), "bb") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 2), "ccc") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 3), "dddd") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 4), "eeeee") == 0);
  for (idx = 0; idx < gt_splitter_size(s) - 1; idx++)
    RESET_DELIMITER(gt_splitter_get_token(s, idx), delimiter2);
  gt_splitter_reset(s);
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split_non_empty(s, string_2, strlen(string_2), delimiter2);
  gt_ensure(gt_splitter_size(s) == 5);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "a") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 1), "bb") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 2), "ccc") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 3), "dddd") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 4), "eeeee") == 0);
  gt_splitter_reset(s);

  /* string_3 */
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split(s, string_3, strlen(string_3), delimiter2);
  gt_ensure(gt_splitter_size(s) == 1);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "") == 0);
  for (idx = 0; idx < gt_splitter_size(s) - 1; idx++)
    RESET_DELIMITER(gt_splitter_get_token(s, idx), delimiter2);
  gt_splitter_reset(s);
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split_non_empty(s, string_3, strlen(string_3), delimiter2);
  gt_ensure(gt_splitter_size(s) == 0);
  gt_splitter_reset(s);

  /* string_4 */
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split(s, string_4, strlen(string_4), delimiter1);
  gt_ensure(gt_splitter_size(s) == 3);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "a") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 1), "") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 2), "b") == 0);
  for (idx = 0; idx < gt_splitter_size(s) - 1; idx++)
    RESET_DELIMITER(gt_splitter_get_token(s, idx), delimiter1);
  gt_splitter_reset(s);
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split_non_empty(s, string_4, strlen(string_4), delimiter1);
  gt_ensure(gt_splitter_size(s) == 2);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "a") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 1), "b") == 0);
  gt_splitter_reset(s);

  /* string_5 */
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split(s, string_5, strlen(string_5), delimiter1);
  gt_ensure(gt_splitter_size(s) == 3);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "ac") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 1), "bc") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 2), "") == 0);
  for (idx = 0; idx < gt_splitter_size(s) - 1; idx++)
    RESET_DELIMITER(gt_splitter_get_token(s, idx), delimiter1);
  gt_splitter_reset(s);
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split_non_empty(s, string_5, strlen(string_5), delimiter1);
  gt_ensure(gt_splitter_size(s) == 2);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "ac") == 0);
  gt_ensure(strcmp(gt_splitter_get_token(s, 1), "bc") == 0);
  gt_splitter_reset(s);

  /* string_6 */
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split(s, string_6, strlen(string_6), delimiter3);
  gt_ensure(gt_splitter_size(s) == 1);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "test") == 0);
  gt_splitter_reset(s);
  gt_ensure(!gt_splitter_size(s));
  gt_splitter_split_non_empty(s, string_6, strlen(string_6), delimiter3);
  gt_ensure(gt_splitter_size(s) == 1);
  gt_ensure(strcmp(gt_splitter_get_token(s, 0), "test") == 0);

  /* free */
  gt_splitter_delete(s);

  return had_err;
}
