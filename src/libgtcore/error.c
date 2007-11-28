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

void error_set(Error *e, const char *format, ...)
{
  va_list ap;
  if (!e) return;
  va_start(ap, format);
  error_vset(e, format, ap);
  va_end(ap);
}

void error_vset(Error *e, const char *format, va_list ap)
{
  assert(e && format);
  e->error_is_set = true;
  (void) vsnprintf(e->error_string, sizeof (e->error_string), format, ap);
}

bool error_is_set(const Error *e)
{
  assert(e);
  return e->error_is_set;
}

void error_unset(Error *e)
{
  assert(e);
  e->error_is_set = false;
  e->error_string[0] = '\0';
}

const char* error_get(const Error *e)
{
  assert(e && e->error_is_set);
  return e->error_string;
}

void error_set_progname(Error *e, const char *progname)
{
  assert(e && progname);
  ma_free(e->progname);
  e->progname = cstr_dup(progname);
}

const char* error_get_progname(const Error *e)
{
  assert(e);
  return e->progname;
}

void error_delete(Error *e)
{
  if (!e) return;
  ma_free(e->progname);
  ma_free(e);
}
