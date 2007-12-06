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

#include <sys/types.h>
#include <assert.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include "libgtcore/ensure.h"
#include "libgtcore/error.h"
#include "libgtcore/grep.h"

static void grep_error(int errcode, regex_t *matcher, Error *err)
{
  char sbuf[BUFSIZ], *buf;
  size_t bufsize;
  error_check(err);
  bufsize = regerror(errcode, matcher, NULL, 0);
  buf = malloc(bufsize);
  (void) regerror(errcode, matcher, buf ? buf : sbuf, buf ? bufsize : BUFSIZ);
  error_set(err, "grep(): %s", buf ? buf : sbuf);
  free(buf);
}

int grep(bool *match, const char *pattern, const char *line, Error *err)
{
  regex_t matcher;
  int rval, had_err = 0;
  error_check(err);
  assert(pattern && line);
  if ((rval = regcomp(&matcher, pattern, REG_EXTENDED | REG_NOSUB))) {
    grep_error(rval, &matcher, err);
    had_err = -1;
  }
  if (!had_err) {
    rval = regexec(&matcher, line, 0, NULL, 0);
    if (rval && rval != REG_NOMATCH) {
      grep_error(rval, &matcher, err);
      had_err = -1;
    }
  }
  regfree(&matcher);
  if (!had_err) {
    if (rval)
      *match = false;
    else
      *match = true;
  }
  return had_err;
}

int grep_unit_test(Error *err)
{
  bool match;
  int grep_err, had_err = 0;
  error_check(err);

  grep_err = grep(&match, "a", "a", NULL);
  ensure(had_err, !grep_err);
  ensure(had_err, match);

  grep_err = grep(&match, "b", "a", NULL);
  ensure(had_err, !grep_err);
  ensure(had_err, !match);

  grep_err =  grep(&match, "aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA", NULL);
  ensure(had_err, !grep_err);
  ensure(had_err, match);

  grep_err = grep(&match, "aba", "wenbapzbpqSayhzzaBaZZqyghaAAahhaA", NULL);
  ensure(had_err, !grep_err);
  ensure(had_err, !match);

  grep_err = grep(&match, "^aba", "abawenbapzbpqSayhzzZZqyghaAAahhaA", NULL);
  ensure(had_err, !grep_err);
  ensure(had_err, match);

  grep_err = grep(&match, "^aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA", NULL);
  ensure(had_err, !grep_err);
  ensure(had_err, !match);

  return had_err;
}
