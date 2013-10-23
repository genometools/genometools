/*
  Copyright (c) 2013 Gordon Gremme <gordon@gremme.org>

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

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "core/compat.h"

int gt_mkstemp(char *templ)
{
#ifndef _WIN32
  return mkstemp(templ);
#else
  /* XXX: is this replacement good enough? */
#ifdef _mktemp_s
  errno_t err = _mktemp_s(templ, strlen(templ) + 1);
  if (err == EINVAL)
    return -1;
#else
  if (!_mktemp(templ))
    return -1;
#endif
  return open(templ, O_RDWR, O_EXCL);
#endif
}

GtUword gt_pagesize(void)
{
#ifndef _WIN32
  return sysconf(_SC_PAGESIZE);
#else
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwPageSize;
#endif
}
