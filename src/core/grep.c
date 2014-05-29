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
#include <stdlib.h>
#include <string.h>
#include "core/ensure.h"
#include "core/grep_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "tre.h"

static void grep_error(int errcode, regex_t *matcher, GtError *err)
{
  char sbuf[BUFSIZ], *buf;
  size_t bufsize;
  gt_error_check(err);
  bufsize = tre_regerror(errcode, matcher, NULL, 0);
  buf = gt_malloc(bufsize);
  (void) tre_regerror(errcode, matcher, buf ? buf : sbuf,
                      buf ? bufsize : BUFSIZ);
  gt_error_set(err, "grep(): %s", buf ? buf : sbuf);
  gt_free(buf);
}

int gt_grep_nt(GT_UNUSED bool *match, GT_UNUSED const char *pattern,
               GT_UNUSED const char *line, size_t len, GtError *err)
{
  regex_t matcher;
  int rval, had_err = 0;
  gt_error_check(err);
  gt_assert(pattern && line);
  if ((rval = tre_regcomp(&matcher, pattern, REG_EXTENDED | REG_NOSUB))) {
    grep_error(rval, &matcher, err);
    had_err = -1;
  }
  if (!had_err) {
    rval = tre_regnexec(&matcher, line, len, 0, NULL, 0);
    if (rval && rval != REG_NOMATCH) {
      grep_error(rval, &matcher, err);
      had_err = -1;
    }
  }
  tre_regfree(&matcher);
  if (!had_err) {
    if (rval)
      *match = false;
    else
      *match = true;
  }
  return had_err;
}

int gt_grep(GT_UNUSED bool *match, GT_UNUSED const char *pattern,
            GT_UNUSED const char *line, GtError *err)
{
  return gt_grep_nt(match, pattern, line, strlen(line), err);
}

int gt_grep_unit_test(GtError *err)
{
  bool match;
  int grep_err, had_err = 0;
  gt_error_check(err);

  grep_err = gt_grep(&match, "a", "a", NULL);
  gt_ensure(!grep_err);
  gt_ensure(match);
  grep_err = gt_grep_nt(&match, "a", "abb", 1, NULL);
  gt_ensure(!grep_err);
  gt_ensure(match);
  grep_err = gt_grep_nt(&match, "b", "abb", 2, NULL);
  gt_ensure(!grep_err);
  gt_ensure(match);

  grep_err = gt_grep(&match, "b", "a", NULL);
  gt_ensure(!grep_err);
  gt_ensure(!match);
  grep_err = gt_grep_nt(&match, "b", "ab", 1, NULL);
  gt_ensure(!grep_err);
  gt_ensure(!match);
  grep_err = gt_grep_nt(&match, "b", "b", 0, NULL);
  gt_ensure(!grep_err);
  gt_ensure(!match);

  grep_err =  gt_grep(&match, "^foo ", "foo bar", NULL);
  gt_ensure(!grep_err);
  gt_ensure(match);
  grep_err =  gt_grep(&match, "^foo ", "baz foo bar", NULL);
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
