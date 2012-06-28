/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2006-2012 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr_api.h"
#include "core/cstr_array.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/xansi_api.h"

char** gt_cstr_split(const char *cstr, char sep)
{
  size_t n = 0;
  char **res;
  unsigned long i,
                lastpos = 0,
                respos = 0;
  gt_assert(cstr);

  for (i = 0; i < strlen(cstr); i++) {
    if (cstr[i] == sep)
      n++;
  }
  n++;
  res = gt_calloc(n+1, sizeof (char*));
  for (i = 0; i < strlen(cstr)+1; i++) {
    if (cstr[i] == sep || cstr[i] == '\0') {
      if (i != 0)
        res[respos] = gt_calloc(i - lastpos + 1, sizeof (char));
      if (i - lastpos >= 1)
        memcpy(res[respos], cstr + lastpos, (i - lastpos) * sizeof (char));
      lastpos = i+1;
      respos++;
    }
  }
  gt_assert(respos == n);
  res[respos] = NULL;
  return res;
}

char* gt_cstr_dup(const char *cstr)
{
  size_t size;
  char *copy;
  gt_assert(cstr);
  size = strlen(cstr) + 1;
  copy = gt_malloc(size);
  memcpy(copy, cstr, size);
  return copy;
}

char* gt_cstr_dup_nt(const char *cstr, unsigned long length)
{
  char *copy;
  gt_assert(cstr);
  copy = gt_malloc(length+1);
  memcpy(copy, cstr, length);
  copy[length] = '\0';
  return copy;
}

void gt_cstr_rep(char *cstr, char f, char t)
{
  char *cc;
  gt_assert(cstr);
  cc = cstr;
  while (*cc) {
    if (*cc == f)
      *cc = t;
    cc++;
  }
}

void gt_cstr_show(const char *cstr, unsigned long length, FILE *fp)
{
  unsigned long i;
  gt_assert(cstr && fp);
  for (i = 0; i < length; i++)
    gt_xfputc(cstr[i], fp);
}

unsigned long gt_cstr_length_up_to_char(const char *cstr, char c)
{
  char *suffix;
  gt_assert(cstr);
  suffix = strchr(cstr, c);
  if (suffix)
    return suffix - cstr;
  return strlen(cstr);
}

char* gt_cstr_rtrim(char *cstr, char remove)
{
  char *tmp;
  gt_assert(cstr);
  for (tmp = cstr + strlen(cstr) - 1; tmp >= cstr && *tmp == remove; --tmp);
  *(tmp+1) = '\0';
  return cstr;
}

int gt_cstr_unit_test(GtError *err)
{
  int had_err = 0;
  char **res;
  const char *foo = "foo  bar baz";
  gt_error_check(err);

  res = gt_cstr_split(foo, ' ');
  gt_ensure(had_err, strcmp(res[0], "foo") == 0);
  gt_ensure(had_err, strcmp(res[1], "") == 0);
  gt_ensure(had_err, strcmp(res[2], "bar") == 0);
  gt_ensure(had_err, strcmp(res[3], "baz") == 0);
  gt_ensure(had_err, res[4] == NULL);
  gt_cstr_array_delete(res);

  res = gt_cstr_split("", ' ');
  gt_ensure(had_err, res[0] == NULL);
  gt_cstr_array_delete(res);

  return had_err;
}
