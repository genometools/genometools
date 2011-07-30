/*
  Copyright (c) 2007-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "core/ma.h"
#include "core/str_array.h"

struct GtStrArray {
  GtArray *strings;
  unsigned int reference_count;
};

GtStrArray* gt_str_array_new(void)
{
  GtStrArray *sa = gt_malloc(sizeof *sa);
  sa->strings = gt_array_new(sizeof (GtStr*));
  sa->reference_count = 0;
  return sa;
}

GtStrArray* gt_str_array_ref(GtStrArray *sa)
{
  if (!sa) return NULL;
  sa->reference_count++;
  return sa;
}

GtStrArray* gt_str_array_new_file(const char *path)
{
  GtStrArray *filecontent;
  GtFile *fpin;
  GtStr *line;
  fpin = gt_file_xopen(path, "r");
  gt_assert(fpin);
  line = gt_str_new();
  filecontent = gt_str_array_new();
  while (gt_str_read_next_line_generic(line, fpin) != EOF) {
    gt_str_array_add_cstr(filecontent, gt_str_get(line));
    gt_str_reset(line);
  }
  gt_str_delete(line);
  gt_file_delete(fpin);
  return filecontent;
}

void gt_str_array_add_cstr(GtStrArray *sa, const char *cstr)
{
  GtStr *str;
  gt_assert(sa && cstr);
  str = gt_str_new_cstr(cstr);
  gt_array_add(sa->strings, str);
}

void gt_str_array_add_cstr_nt(GtStrArray *sa, const char *cstr,
                             unsigned long length)
{
  GtStr *str;
  gt_assert(sa && cstr);
  str = gt_str_new();
  gt_str_append_cstr_nt(str, cstr, length);
  gt_array_add(sa->strings, str);
}

void gt_str_array_add(GtStrArray *sa, const GtStr *str)
{
  GtStr *clone;
  gt_assert(sa && str);
  clone = gt_str_clone(str);
  gt_array_add(sa->strings, clone);
}

const char* gt_str_array_get(const GtStrArray *sa, unsigned long strnum)
{
  gt_assert(sa && strnum < gt_array_size(sa->strings));
  return gt_str_get(*(GtStr**) gt_array_get(sa->strings, strnum));
}

GtStr* gt_str_array_get_str(const GtStrArray *sa, unsigned long strnum)
{
  gt_assert(sa && strnum < gt_array_size(sa->strings));
  return *(GtStr**) gt_array_get(sa->strings, strnum);
}

void gt_str_array_set_cstr(GtStrArray *sa, unsigned long strnum,
                           const char *cstr)
{
  GtStr *str;
  gt_assert(sa && strnum < gt_array_size(sa->strings) && cstr);
  str = *(GtStr**) gt_array_get(sa->strings, strnum);
  gt_str_set(str, cstr);
}

void gt_str_array_set(GtStrArray *sa, unsigned long strnum, const GtStr *instr)
{
  GtStr *str;
  gt_assert(sa && strnum < gt_array_size(sa->strings) && instr);
  str = *(GtStr**) gt_array_get(sa->strings, strnum);
  gt_str_set(str, gt_str_get(instr));
}

void gt_str_array_set_size(GtStrArray *sa, unsigned long size)
{
  unsigned long i;
  gt_assert(sa && size <= gt_array_size(sa->strings));
  for (i = size; i < gt_array_size(sa->strings); i++)
    gt_str_delete(*(GtStr**) gt_array_get(sa->strings, i));
  gt_array_set_size(sa->strings, size);
}

void gt_str_array_reset(GtStrArray *sa)
{
  gt_assert(sa);
  gt_str_array_set_size(sa, 0);
}

unsigned long gt_str_array_size(const GtStrArray *sa)
{
  gt_assert(sa);
  return gt_array_size(sa->strings);
}

void gt_str_array_delete(GtStrArray *sa)
{
  unsigned long i;
  if (!sa) return;
  if (sa->reference_count) {
    sa->reference_count--;
    return;
  }
  for (i = 0; i < gt_array_size(sa->strings); i++)
    gt_str_delete(*(GtStr**) gt_array_get(sa->strings, i));
  gt_array_delete(sa->strings);
  gt_free(sa);
}
