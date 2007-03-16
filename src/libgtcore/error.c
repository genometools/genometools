/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdarg.h>
#include <libgtcore/env.h>
#include <libgtcore/xansi.h>

struct Error {
  char error_string[BUFSIZ];
  bool error_is_set;
};

Error* error_new(MA *ma)
{
  return ma_calloc(ma, 1, sizeof (Error));
}

void error_set(Error *e, const char *format, ...)
{
  va_list ap;
  assert(e);
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

void error_delete(Error *e, MA *ma)
{
  if (!e) return;
  ma_free(e, ma);
}
