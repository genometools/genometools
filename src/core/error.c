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
#include "core/cstr_api.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/xansi_api.h"

struct GtError {
  char error_string[BUFSIZ],
       *progname;
  bool error_is_set;
};

GtError* gt_error_new(void)
{
  return gt_calloc(1, sizeof (GtError));
}

void gt_error_set(GtError *err, const char *format, ...)
{
  va_list ap;
  if (!err) return;
  va_start(ap, format);
  gt_error_vset(err, format, ap);
  va_end(ap);
}

void gt_error_set_nonvariadic(GtError *err, const char *msg)
{
  if (!err) return;
  gt_error_set(err, "%s", msg);
}

void gt_error_vset(GtError *err, const char *format, va_list ap)
{
  gt_assert(err && format);
  err->error_is_set = true;
  (void) vsnprintf(err->error_string, sizeof (err->error_string), format, ap);
}

bool gt_error_is_set(const GtError *err)
{
  gt_assert(err);
  return err->error_is_set;
}

void gt_error_unset(GtError *err)
{
  gt_assert(err);
  err->error_is_set = false;
  err->error_string[0] = '\0';
}

const char* gt_error_get(const GtError *err)
{
  gt_assert(err && err->error_is_set);
  return err->error_string;
}

void gt_error_set_progname(GtError *err, const char *progname)
{
  gt_assert(err && progname);
  gt_free(err->progname);
  err->progname = gt_cstr_dup(progname);
}

const char* gt_error_get_progname(const GtError *err)
{
  gt_assert(err);
  return err->progname;
}

void gt_error_delete(GtError *err)
{
  if (!err) return;
  gt_free(err->progname);
  gt_free(err);
}
