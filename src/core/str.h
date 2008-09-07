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

#ifndef STR_H
#define STR_H

#include <stdio.h>
#include "core/error.h"
#include "core/genfile.h"

/* the string class, string objects are strings which grow on demand */
typedef struct GT_Str GT_Str;

GT_Str*       str_new(void);
GT_Str*       str_new_cstr(const char*);
GT_Str*       str_clone(const GT_Str*);
GT_Str*       str_ref(GT_Str*);
/* never returns NULL, always '\0' terminated */
char*         str_get(const GT_Str*);
/* never returns NULL, not always '\0' terminated */
void*         str_get_mem(const GT_Str*);
void          str_set(GT_Str*, const char*);
void          str_append_str(GT_Str*, const GT_Str*);
void          str_append_cstr(GT_Str*, const char*);
/* appends an unterminated cstr */
void          str_append_cstr_nt(GT_Str*, const char*, unsigned long);
void          str_append_char(GT_Str*, char);
void          str_append_double(GT_Str*, double, int precision);
void          str_append_ulong(GT_Str*, unsigned long);
void          str_set_length(GT_Str*, unsigned long);
void          str_reset(GT_Str*);
int           str_cmp(const GT_Str*, const GT_Str*);
int           str_read_next_line(GT_Str*, FILE*);
int           str_read_next_line_generic(GT_Str*, GenFile*);
unsigned long str_length(const GT_Str*);
int           str_unit_test(GT_Error*);
void          str_delete(GT_Str*);

#endif
