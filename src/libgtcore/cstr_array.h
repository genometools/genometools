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

#ifndef CSTR_ARRAY_H
#define CSTR_ARRAY_H

#include "libgtcore/genfile.h"

/* Return copy of <cstr_array>. */
char**        cstr_array_dup(const char **cstr_array);
/* Use p and a blank as prefix for cstr_array[0] and return the result. */
char**        cstr_array_prefix_first(const char **cstr_array, const char *p);
char**        cstr_array_preprend(const char **cstr_array, const char *p);
void          cstr_array_show(char **cstr_array, FILE*);
void          cstr_array_show_genfile(const char **cstr_array, GenFile*);
unsigned long cstr_array_size(const char **cstr_array); /* O(n) */
void          cstr_array_delete(char **cstr_array);

#endif
