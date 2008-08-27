/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef CSTR_TABLE_H
#define CSTR_TABLE_H

#include "libgtcore/error.h"

/* A table of C-strings. */
typedef struct CstrTable CstrTable;

CstrTable*  cstr_table_new();
void        cstr_table_delete(CstrTable*);
/* Add <cstr> to <table>. */
void        cstr_table_add(CstrTable *table, const char *cstr);
/* If a C-string equal to <cstr> is contained in <table>, it is returned.
   Otherwise NULL is returned. */
const char* cstr_table_get(const CstrTable *table, const char *cstr);
int         cstr_table_unit_test(Error*);

#endif
