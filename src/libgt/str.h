/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STR_H
#define STR_H

#include <stdio.h>
#include "env.h"

/* the string class, string objects are strings which grow on demand */
typedef struct Str Str;

Str*          str_new(Env*);
Str*          str_new_cstr(const char*, Env*);
Str*          str_clone(const Str*, Env*);
Str*          str_ref(Str*);
char*         str_get(const Str*); /* never returns NULL */
void          str_set(Str*, const char*, Env*);
void          str_append_str(Str*, const Str*, Env*);
void          str_append_cstr(Str*, const char*, Env*);
/* appends an unterminated cstr */
void          str_append_cstr_nt(Str*, const char*, unsigned long, Env*);
void          str_append_char(Str*, char, Env*);
void          str_append_ulong(Str*, unsigned long, Env*);
void          str_set_length(Str*, unsigned long);
void          str_reset(Str*);
int           str_cmp(const Str*, const Str*);
int           str_read_next_line(Str*, FILE*, Env*);
unsigned long str_length(const Str*);
int           str_unit_test(Env*);
void          str_delete(Str*, Env*);

#endif
