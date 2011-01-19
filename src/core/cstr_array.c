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
#include "core/cstr_array.h"
#include "core/ma.h"
#include "core/xansi_api.h"

char** gt_cstr_array_dup(const char **gt_cstr_array)
{
  unsigned long i, size = 0;
  char **copy;
  gt_assert(gt_cstr_array);
  while (gt_cstr_array[size++]);
  copy = gt_malloc(size * sizeof *copy);
  for (i = 0; i < size - 1; i++)
    copy[i] = gt_cstr_dup(gt_cstr_array[i]);
  copy[size-1] = NULL;
  return copy;
}

char** gt_cstr_array_prefix_first(const char **gt_cstr_array, const char *p)
{
  unsigned long i, a_len, f_len;
  char **a;
  gt_assert(gt_cstr_array && p);
  a_len = gt_cstr_array_size(gt_cstr_array);
  a = gt_malloc(sizeof (char*) * (a_len + 1));
  f_len = strlen(p) + strlen(gt_cstr_array[0]) + 2; /* blank + '\0' */
  a[0] = gt_malloc(sizeof (char) * f_len);
  (void) snprintf(a[0], f_len, "%s %s", p, gt_cstr_array[0]);
  for (i = 1; i < a_len; i++)
    a[i] = gt_cstr_dup(gt_cstr_array[i]);
  a[a_len] = NULL;
  return a;
}

char** gt_cstr_array_preprend(const char **gt_cstr_array, const char *p)
{
  unsigned long i, a_len;
  char **a;
  gt_assert(gt_cstr_array && p);
  a_len = gt_cstr_array_size(gt_cstr_array);
  a = gt_malloc(sizeof (char*) * (a_len + 2));
  a[0] = gt_cstr_dup(p);
  for (i = 0; i < a_len; i++)
    a[i+1] = gt_cstr_dup(gt_cstr_array[i]);
  a[a_len+1] = NULL;
  return a;
}

void gt_cstr_array_show(char **gt_cstr_array, FILE *fp)
{
  unsigned long i = 0;
  while (gt_cstr_array[i]) {
    if (i)
      gt_xfputc(' ', fp);
    gt_xfputs(gt_cstr_array[i], fp);
    i++;
  }
  gt_xfputc('\n', fp);
}

void gt_cstr_array_show_genfile(const char **gt_cstr_array, GtFile *fp)
{
  unsigned long i = 0;
  while (gt_cstr_array[i]) {
    if (i)
      gt_file_xfputc(' ', fp);
    gt_file_xfputs(gt_cstr_array[i], fp);
    i++;
  }
  gt_file_xfputc('\n', fp);
}

unsigned long gt_cstr_array_size(const char **gt_cstr_array)
{
  unsigned long i = 0;
  while (gt_cstr_array[i])
    i++;
  return i;
}

void gt_cstr_array_delete(char **gt_cstr_array)
{
  unsigned long i = 0;
  if (!gt_cstr_array) return;
  while (gt_cstr_array[i])
    gt_free(gt_cstr_array[i++]);
  gt_free(gt_cstr_array);
}
