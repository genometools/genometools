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

#include <string.h>
#include "core/cstr_api.h"
#include "core/ma.h"
#include "core/xansi_api.h"

char* gt_cstr_dup(const char *cstr)
{
  size_t size;
  char *copy;
  gt_assert(cstr);
  size = strlen(cstr) + 1;
  copy = gt_malloc(size);
  memcpy(copy, cstr, size);
  return copy;
}

char* gt_cstr_dup_nt(const char *cstr, unsigned long length)
{
  char *copy;
  gt_assert(cstr);
  copy = gt_malloc(length+1);
  memcpy(copy, cstr, length);
  copy[length] = '\0';
  return copy;
}

void gt_cstr_rep(char *cstr, char f, char t)
{
  char *cc;
  gt_assert(cstr);
  cc = cstr;
  while (*cc) {
    if (*cc == f)
      *cc = t;
    cc++;
  }
}

void gt_cstr_show(const char *cstr, unsigned long length, FILE *fp)
{
  unsigned long i;
  gt_assert(cstr && fp);
  for (i = 0; i < length; i++)
    gt_xfputc(cstr[i], fp);
}

unsigned long gt_cstr_length_up_to_char(const char *cstr, char c)
{
  char *suffix;
  gt_assert(cstr);
  suffix = strchr(cstr, c);
  if (suffix)
    return suffix - cstr;
  return strlen(cstr);
}

char* gt_cstr_rtrim(char *cstr, char remove)
{
  char *tmp;
  gt_assert(cstr);
  for (tmp = cstr + strlen(cstr) - 1; tmp >= cstr && *tmp == remove; --tmp);
  *(tmp+1) = '\0';
  return cstr;
}
