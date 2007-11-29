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

#include <assert.h>
#include <stdlib.h>
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

char** cstr_array_prefix_first(const char **cstr_array, const char *p)
{
  unsigned long i, a_len, f_len;
  char **a;
  assert(cstr_array && p);
  a_len = cstr_array_size(cstr_array);
  a = ma_malloc(sizeof (char*) * (a_len + 1));
  f_len = strlen(p) + strlen(cstr_array[0]) + 2; /* blank + '\0' */
  a[0] = ma_malloc(sizeof (char) * f_len);
  (void) snprintf(a[0], f_len, "%s %s", p, cstr_array[0]);
  for (i = 1; i < a_len; i++)
    a[i] = cstr_dup(cstr_array[i]);
  a[a_len] = NULL;
  return a;
}

char** cstr_array_preprend(const char **cstr_array, const char *p)
{
  unsigned long i, a_len;
  char **a;
  assert(cstr_array && p);
  a_len = cstr_array_size(cstr_array);
  a = ma_malloc(sizeof (char*) * (a_len + 2));
  a[0] = cstr_dup(p);
  for (i = 0; i < a_len; i++)
    a[i+1] = cstr_dup(cstr_array[i]);
  a[a_len+1] = NULL;
  return a;
}

void cstr_array_show(char **cstr_array, FILE *fp)
{
  unsigned long i = 0;
  while (cstr_array[i]) {
    if (i)
      xfputc(' ', fp);
    xfputs(cstr_array[i], fp);
    i++;
  }
  xfputc('\n', fp);
}

void cstr_array_show_genfile(const char **cstr_array, GenFile *fp)
{
  unsigned long i = 0;
  while (cstr_array[i]) {
    if (i)
      genfile_xfputc(' ', fp);
    genfile_xfputs(cstr_array[i], fp);
    i++;
  }
  genfile_xfputc('\n', fp);
}

unsigned long cstr_array_size(const char **cstr_array)
{
  unsigned long i = 0;
  while (cstr_array[i])
    i++;
  return i;
}

void cstr_array_delete(char **cstr_array)
{
  unsigned long i = 0;
  if (!cstr_array) return;
  while (cstr_array[i])
    ma_free(cstr_array[i++]);
  ma_free(cstr_array);
}
