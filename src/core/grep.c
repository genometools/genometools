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

#ifndef S_SPLINT_S
#include <sys/types.h>
#endif
#ifndef _WIN32
#include <regex.h>
#endif
#include <stdlib.h>
#include "core/ensure.h"
#include "core/grep_api.h"
#include "core/unused_api.h"

#ifndef _WIN32
static void grep_error(int errcode, regex_t *matcher, GtError *err)
{
  char sbuf[BUFSIZ], *buf;
  size_t bufsize;
  gt_error_check(err);
  bufsize = regerror(errcode, matcher, NULL, 0);
  buf = malloc(bufsize);
  (void) regerror(errcode, matcher, buf ? buf : sbuf, buf ? bufsize : BUFSIZ);
  gt_error_set(err, "grep(): %s", buf ? buf : sbuf);
  free(buf);
}
#endif

int gt_grep(GT_UNUSED bool *match, GT_UNUSED const char *pattern,
            GT_UNUSED const char *line, GtError *err)
{
#ifndef _WIN32
  regex_t matcher;
  int rval, had_err = 0;
  gt_error_check(err);
  gt_assert(pattern && line);
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
#else
  /* XXX */
  gt_error_set(err, "gt_grep() not implemented");
  return -1;
#endif
}

int gt_grep_unit_test(GtError *err)
{
  bool match;
  int grep_err, had_err = 0;
  gt_error_check(err);

  grep_err = gt_grep(&match, "a", "a", NULL);
  gt_ensure(!grep_err);
  gt_ensure(match);

  grep_err = gt_grep(&match, "b", "a", NULL);
  gt_ensure(!grep_err);
  gt_ensure(!match);

  grep_err =  gt_grep(&match, "aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA", NULL);
  gt_ensure(!grep_err);
  gt_ensure(match);

  grep_err = gt_grep(&match, "aba", "wenbapzbpqSayhzzaBaZZqyghaAAahhaA", NULL);
  gt_ensure(!grep_err);
  gt_ensure(!match);

  grep_err = gt_grep(&match, "^aba", "abawenbapzbpqSayhzzZZqyghaAAahhaA", NULL);
  gt_ensure(!grep_err);
  gt_ensure(match);

  grep_err = gt_grep(&match, "^aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA", NULL);
  gt_ensure(!grep_err);
  gt_ensure(!match);

  return had_err;
}
