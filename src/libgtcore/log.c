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
  fprintf(l->logfp, "debug: ");
  (void) vfprintf(l->logfp, format, ap);
  (void) putc('\n', l->logfp);
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
