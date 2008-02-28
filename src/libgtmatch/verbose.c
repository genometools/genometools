/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <assert.h>
#include <stdarg.h>
#include "libgtcore/ma.h"
#include "verbose-def.h"

struct Verboseinfo
{
  bool beverbose;
};

void showdefinitelyverbose(const char *format, ...)
{
  va_list ap;

  assert(format != NULL);
  va_start(ap, format);
  printf("# ");
  (void) vprintf(format, ap);
  printf("\n");
  va_end(ap);
}

void showverbose(Verboseinfo *verboseinfo,const char *format, ...)
{
  if (verboseinfo != NULL && verboseinfo->beverbose)
  {
    va_list ap;

    assert(format != NULL);
    va_start(ap, format);
    printf("# ");
    (void) vprintf(format, ap);
    (void) fflush(stdout);
    printf("\n");
    va_end(ap);
  }
}

Verboseinfo *newverboseinfo(bool verbose)
{
  Verboseinfo *v;

  v = ma_malloc(sizeof (Verboseinfo));
  v->beverbose = verbose;
  return v;
}

void freeverboseinfo(Verboseinfo **v)
{
  ma_free(*v);
  *v = NULL;
}
