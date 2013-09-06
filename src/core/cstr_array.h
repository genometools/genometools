/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef CSTR_ARRAY_H
#define CSTR_ARRAY_H

#include "core/file.h"
#include "core/types_api.h"

/* Return copy of <gt_cstr_array>. */
char**  gt_cstr_array_dup(const char **cstr_array);
/* Use p and a blank as prefix for cstr_array[0] and return the result. */
char**  gt_cstr_array_prefix_first(const char **cstr_array, const char *p);
char**  gt_cstr_array_preprend(const char **cstr_array, const char *p);
void    gt_cstr_array_show(char **cstr_array, FILE*);
void    gt_cstr_array_show_genfile(const char **cstr_array, GtFile*);
GtUword gt_cstr_array_size(const char **cstr_array); /* O(n) */
void    gt_cstr_array_delete(char **cstr_array);

#endif
