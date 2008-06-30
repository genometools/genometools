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

#include <assert.h>
#include <string.h>
#include "libgtcore/cstr.h"
#include "libgtcore/ma.h"
#include "libgtcore/xansi.h"

char* cstr_dup(const char *cstr)
{
  size_t size;
  char *copy;
  assert(cstr);
  size = strlen(cstr) + 1;
  copy = ma_malloc(size);
  memcpy(copy, cstr, size);
  return copy;
}

void cstr_rep(char *cstr, char f, char t)
{
  char *cc;
  assert(cstr);
  cc = cstr;
  while (*cc) {
    if (*cc == f)
      *cc = t;
    cc++;
  }
}

void cstr_show(const char *cstr, unsigned long length, FILE *fp)
{
  unsigned long i;
  assert(cstr && fp);
  for (i = 0; i < length; i++)
    xfputc(cstr[i], fp);
}

unsigned long cstr_length_up_to_char(const char *cstr, char c)
{
  char *suffix;
  assert(cstr);
  suffix = strchr(cstr, c);
  if (suffix)
    return suffix - cstr;
  return strlen(cstr);
}
