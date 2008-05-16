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

#include "libgtcore/array.h"
#include "libgtcore/cstr.h"
#include "libgtcore/ma.h"
#include "libgtcore/strarray.h"

struct StrArray {
  Array *strings;
};

StrArray* strarray_new(void)
{
  StrArray *sa;
  sa = ma_malloc(sizeof (StrArray));
  sa->strings = array_new(sizeof (Str*));
  return sa;
}

StrArray* strarray_new_file(const char *path)
{
  StrArray *filecontent;
  GenFile *fpin;
  Str *line;
  fpin = genfile_xopen(path, "r");
  assert(fpin);
  line = str_new();
  filecontent = strarray_new();
  while (str_read_next_line_generic(line, fpin) != EOF) {
    strarray_add_cstr(filecontent, str_get(line));
    str_reset(line);
  }
  str_delete(line);
  genfile_close(fpin);
  return filecontent;
}

void strarray_add_cstr(StrArray *sa, const char *cstr)
{
  Str *str;
  assert(sa && cstr);
  str = str_new_cstr(cstr);
  array_add(sa->strings, str);
}

void strarray_add(StrArray *sa, const Str *str)
{
  Str *clone;
  assert(sa && str);
  clone = str_clone(str);
  array_add(sa->strings, clone);
}

const char* strarray_get(const StrArray *sa, unsigned long strnum)
{
  assert(sa && strnum < array_size(sa->strings));
  return str_get(*(Str**) array_get(sa->strings, strnum));
}

Str* strarray_get_str(const StrArray *sa, unsigned long strnum)
{
  assert(sa && strnum < array_size(sa->strings));
  return *(Str**) array_get(sa->strings, strnum);
}

void strarray_set_size(StrArray *sa, unsigned long size)
{
  unsigned long i;
  assert(sa && size <= array_size(sa->strings));
  for (i = size; i < array_size(sa->strings); i++)
    str_delete(*(Str**) array_get(sa->strings, i));
  array_set_size(sa->strings, size);
}

unsigned long strarray_size(const StrArray *sa)
{
  assert(sa);
  return array_size(sa->strings);
}

void strarray_delete(StrArray *sa)
{
  unsigned long i;
  if (!sa) return;
  for (i = 0; i < array_size(sa->strings); i++)
    str_delete(*(Str**) array_get(sa->strings, i));
  array_delete(sa->strings);
  ma_free(sa);
}
