/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/array.h"
#include "libgtcore/cstr.h"
#include "libgtcore/ma.h"
#include "libgtcore/strarray.h"

struct StrArray {
  Array *strings;
};

StrArray* strarray_new(Env *env)
{
  StrArray *sa;
  env_error_check(env);
  sa = ma_malloc(sizeof (StrArray));
  sa->strings = array_new(sizeof (Str*), env);
  return sa;
}

StrArray* strarray_new_file(const char *path, Env *env)
{
  StrArray *filecontent;
  GenFile *fpin;
  Str *line;
  env_error_check(env);
  fpin = genfile_xopen(path, "r", env);
  assert(fpin);
  line = str_new();
  filecontent = strarray_new(env);
  while (str_read_next_line_generic(line, fpin, env) != EOF) {
    strarray_add_cstr(filecontent,str_get(line),env);
    str_reset(line);
  }
  str_delete(line);
  genfile_xclose(fpin, env);
  return filecontent;
}

void strarray_add_cstr(StrArray *sa, const char *cstr, Env *env)
{
  Str *str;
  env_error_check(env);
  assert(sa && cstr);
  str = str_new_cstr(cstr, env);
  array_add(sa->strings, str, env);
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

unsigned long strarray_size(const StrArray *sa)
{
  assert(sa);
  return array_size(sa->strings);
}

void strarray_delete(StrArray *sa, Env *env)
{
  unsigned long i;
  if (!sa) return;
  for (i = 0; i < array_size(sa->strings); i++)
    str_delete(*(Str**) array_get(sa->strings, i));
  array_delete(sa->strings, env);
  ma_free(sa);
}
