/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdarg.h>
#include <stdio.h>
#include <libgtcore/env.h>
#include <libgtcore/xansi.h>

struct Log {
  FILE *logfp;
};

Log* log_new(MA *ma)
{
  Log *l = ma_malloc(ma, sizeof (Log));
  l->logfp = stderr;
  return l;
}

void log_log(Log *l, const char *format, ...)
{
  va_list ap;
  if (!l) return;
  va_start(ap, format);
  log_vlog(l, format, ap);
  va_end(ap);
}

void log_vlog(Log *l, const char *format, va_list ap)
{
  if (!l) return;
  fprintf(stderr, "debug: ");
  (void) vfprintf(stderr, format, ap);
  (void) putc('\n', stderr);
}

FILE* log_fp(Log *l)
{
  assert(l);
  return l->logfp;
}

void log_delete(Log *l, MA *ma)
{
  if (!l) return;
  ma_free(l, ma);
}
