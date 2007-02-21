/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdarg.h>
#include <stdio.h>
#include "log.h"
#include "xansi.h"

struct Log {
  FILE *logfp;
};

Log* log_new(void)
{
  Log *l = xmalloc(sizeof (Log));
  l->logfp = stderr;
  return l;
}

void log_log(Log *l, const char *format, ...)
{
  va_list ap;
  if (!l) return;
  va_start(ap, format);
  fprintf(stderr, "debug: ");
  (void) vfprintf(stderr, format, ap);
  (void) putc('\n', stderr);
  va_end(ap);
}

FILE* log_fp(Log *l)
{
  assert(l);
  return l->logfp;
}

void log_delete(Log *l)
{
  if (!l) return;
  free(l);
}
