/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/array.h"
#include "core/cstr.h"
#include "core/ma.h"
#include "core/strarray.h"

struct GT_StrArray {
  GT_Array *strings;
};

GT_StrArray* gt_strarray_new(void)
{
  GT_StrArray *sa = ma_malloc(sizeof *sa);
  sa->strings = gt_array_new(sizeof (Str*));
  return sa;
}

GT_StrArray* gt_strarray_new_file(const char *path)
{
  GT_StrArray *filecontent;
  GenFile *fpin;
  Str *line;
  fpin = genfile_xopen(path, "r");
  assert(fpin);
  line = str_new();
  filecontent = gt_strarray_new();
  while (str_read_next_line_generic(line, fpin) != EOF) {
    gt_strarray_add_cstr(filecontent, str_get(line));
    str_reset(line);
  }
  str_delete(line);
  genfile_close(fpin);
  return filecontent;
}

void gt_strarray_add_cstr(GT_StrArray *sa, const char *cstr)
{
  Str *str;
  assert(sa && cstr);
  str = str_new_cstr(cstr);
  gt_array_add(sa->strings, str);
}

void gt_strarray_add_cstr_nt(GT_StrArray *sa, const char *cstr,
                             unsigned long length)
{
  Str *str;
  assert(sa && cstr);
  str = str_new();
  str_append_cstr_nt(str, cstr, length);
  gt_array_add(sa->strings, str);
}

void gt_strarray_add(GT_StrArray *sa, const Str *str)
{
  Str *clone;
  assert(sa && str);
  clone = str_clone(str);
  gt_array_add(sa->strings, clone);
}

const char* gt_strarray_get(const GT_StrArray *sa, unsigned long strnum)
{
  assert(sa && strnum < gt_array_size(sa->strings));
  return str_get(*(Str**) gt_array_get(sa->strings, strnum));
}

Str* gt_strarray_get_str(const GT_StrArray *sa, unsigned long strnum)
{
  assert(sa && strnum < gt_array_size(sa->strings));
  return *(Str**) gt_array_get(sa->strings, strnum);
}

void gt_strarray_set_size(GT_StrArray *sa, unsigned long size)
{
  unsigned long i;
  assert(sa && size <= gt_array_size(sa->strings));
  for (i = size; i < gt_array_size(sa->strings); i++)
    str_delete(*(Str**) gt_array_get(sa->strings, i));
  gt_array_set_size(sa->strings, size);
}

unsigned long gt_strarray_size(const GT_StrArray *sa)
{
  assert(sa);
  return gt_array_size(sa->strings);
}

void gt_strarray_delete(GT_StrArray *sa)
{
  unsigned long i;
  if (!sa) return;
  for (i = 0; i < gt_array_size(sa->strings); i++)
    str_delete(*(Str**) gt_array_get(sa->strings, i));
  gt_array_delete(sa->strings);
  ma_free(sa);
}
