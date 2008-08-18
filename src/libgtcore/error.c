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
#include "libgtcore/cstr.h"
#include "libgtcore/error.h"
#include "libgtcore/ma.h"
#include "libgtcore/xansi.h"

struct Error {
  char error_string[BUFSIZ],
       *progname;
  bool error_is_set;
};

Error* error_new(void)
{
  return ma_calloc(1, sizeof (Error));
}

void error_set(Error *err, const char *format, ...)
{
  va_list ap;
  if (!err) return;
  va_start(ap, format);
  error_vset(err, format, ap);
  va_end(ap);
}

void error_vset(Error *err, const char *format, va_list ap)
{
  assert(err && format);
  err->error_is_set = true;
  (void) vsnprintf(err->error_string, sizeof (err->error_string), format, ap);
}

bool error_is_set(const Error *err)
{
  assert(err);
  return err->error_is_set;
}

void error_unset(Error *err)
{
  assert(err);
  err->error_is_set = false;
  err->error_string[0] = '\0';
}

const char* error_get(const Error *err)
{
  assert(err && err->error_is_set);
  return err->error_string;
}

void error_set_progname(Error *err, const char *progname)
{
  assert(err && progname);
  ma_free(err->progname);
  err->progname = cstr_dup(progname);
}

const char* error_get_progname(const Error *err)
{
  assert(err);
  return err->progname;
}

void error_delete(Error *err)
{
  if (!err) return;
  ma_free(err->progname);
  ma_free(err);
}
