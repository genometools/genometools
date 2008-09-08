/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <stdarg.h>
#include "core/cstr.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/xansi.h"

struct GT_Error {
  char gt_error_string[BUFSIZ],
       *progname;
  bool gt_error_is_set;
};

GT_Error* gt_error_new(void)
{
  return gt_calloc(1, sizeof (GT_Error));
}

void gt_error_set(GT_Error *err, const char *format, ...)
{
  va_list ap;
  if (!err) return;
  va_start(ap, format);
  gt_error_vset(err, format, ap);
  va_end(ap);
}

void gt_error_vset(GT_Error *err, const char *format, va_list ap)
{
  assert(err && format);
  err->gt_error_is_set = true;
  (void) vsnprintf(err->gt_error_string, sizeof (err->gt_error_string), format, ap);
}

bool gt_error_is_set(const GT_Error *err)
{
  assert(err);
  return err->gt_error_is_set;
}

void gt_error_unset(GT_Error *err)
{
  assert(err);
  err->gt_error_is_set = false;
  err->gt_error_string[0] = '\0';
}

const char* gt_error_get(const GT_Error *err)
{
  assert(err && err->gt_error_is_set);
  return err->gt_error_string;
}

void gt_error_set_progname(GT_Error *err, const char *progname)
{
  assert(err && progname);
  gt_free(err->progname);
  err->progname = cstr_dup(progname);
}

const char* gt_error_get_progname(const GT_Error *err)
{
  assert(err);
  return err->progname;
}

void gt_error_delete(GT_Error *err)
{
  if (!err) return;
  gt_free(err->progname);
  gt_free(err);
}
