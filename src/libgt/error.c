/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdarg.h>
#include "error.h"
#include "xansi.h"

/* the simple interface */
void error(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);
  fprintf(stderr, "error: ");
  (void) vfprintf(stderr, format, ap);
  (void) putc('\n', stderr);
  va_end(ap);
  exit (EXIT_FAILURE);
}

/* the sophisticated interface */
struct Error {
  char error_string[BUFSIZ];
  bool error_is_set;
};

Error* error_new(void)
{
  return xcalloc(1, sizeof(Error));
}

void error_set(Error *e, const char *format, ...)
{
  va_list ap;
  assert(e);
  va_start(ap, format);
  e->error_is_set = true;
  (void) vsnprintf(e->error_string, sizeof(e->error_string), format, ap);
  va_end(ap);
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

void error_abort(const Error *e)
{
  assert(e);
  if (e->error_is_set) {
    fprintf(stderr, "error: %s\n", e->error_string);
    exit(EXIT_FAILURE);
  }
}

void error_free(Error *e)
{
  if (!e) return;
  free(e);
}
