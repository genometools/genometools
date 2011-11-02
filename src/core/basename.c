/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "core/ensure.h"
#include "core/basename_api.h"
#include "core/ma.h"

char *gt_basename(const char *path)
{
  char *sbuf, *c;
  bool foundother = false;
  size_t i, pathlen;

  if (path != NULL)
    pathlen = strlen(path);
  else
    pathlen = 0;
  sbuf = gt_malloc(sizeof (char) * (pathlen + 2));
  if (path == NULL || *path == '\0') {
    strcpy(sbuf, ".");
    return sbuf;
  }
  strcpy(sbuf, path);
  for (c = sbuf + pathlen - 1; c >= sbuf; c--) {
    if (*c == '/') {
      if (foundother) {
        c++;
        for (i=0; c[i] != '\0'; i++)
          sbuf[i] = c[i];
        sbuf[i] = '\0';
        break;
      }
      if (c > sbuf)
        *c = '\0';
    }
    else
      foundother = true;
  }
  return sbuf;
}

int gt_basename_unit_test(GtError *err)
{
  char *bn;
  int had_err = 0;
  gt_error_check(err);

  bn = gt_basename("/usr/lib");
  gt_ensure(had_err, !strcmp(bn, "lib"));
  gt_free(bn);

  if (!had_err) {
    bn = gt_basename("/usr/");
    gt_ensure(had_err, !strcmp(bn, "usr"));
    gt_free(bn);
  }

  if (!had_err) {
    bn = gt_basename("usr");
    gt_ensure(had_err, !strcmp(bn, "usr"));
    gt_free(bn);
  }

  if (!had_err) {
    bn = gt_basename("/");
    gt_ensure(had_err, !strcmp(bn, "/"));
    gt_free(bn);
  }

  if (!had_err) {
    bn = gt_basename("///");
    gt_ensure(had_err, !strcmp(bn, "/"));
    gt_free(bn);
  }

  if (!had_err) {
    bn = gt_basename("//usr//lib//");
    gt_ensure(had_err, !strcmp(bn, "lib"));
    gt_free(bn);
  }

  if (!had_err) {
    bn = gt_basename(NULL);
    gt_ensure(had_err, !strcmp(bn, "."));
    gt_free(bn);
  }

  if (!had_err) {
    bn = gt_basename("");
    gt_ensure(had_err, !strcmp(bn, "."));
    gt_free(bn);
  }

  if (!had_err) {
    bn = gt_basename(".");
    gt_ensure(had_err, !strcmp(bn, "."));
    gt_free(bn);
  }

  if (!had_err) {
    bn = gt_basename("..");
    gt_ensure(had_err, !strcmp(bn, ".."));
    gt_free(bn);
  }

  return had_err;
}
