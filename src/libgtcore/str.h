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
#include "libgtcore/error.h"
#include "libgtcore/genfile.h"

/* the string class, string objects are strings which grow on demand */
typedef struct Str Str;

Str*          str_new(void);
Str*          str_new_cstr(const char*);
Str*          str_clone(const Str*);
Str*          str_ref(Str*);
/* never returns NULL, always '\0' terminated */
char*         str_get(const Str*);
/* never returns NULL, not always '\0' terminated */
void*         str_get_mem(const Str*);
void          str_set(Str*, const char*);
void          str_append_str(Str*, const Str*);
void          str_append_cstr(Str*, const char*);
/* appends an unterminated cstr */
void          str_append_cstr_nt(Str*, const char*, unsigned long);
void          str_append_char(Str*, char);
void          str_append_double(Str*, double, int precision);
void          str_append_ulong(Str*, unsigned long);
void          str_set_length(Str*, unsigned long);
void          str_reset(Str*);
int           str_cmp(const Str*, const Str*);
int           str_read_next_line(Str*, FILE*);
int           str_read_next_line_generic(Str*, GenFile*);
unsigned long str_length(const Str*);
int           str_unit_test(Error*);
void          str_delete(Str*);

#endif
