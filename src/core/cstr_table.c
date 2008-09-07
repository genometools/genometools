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

#include "core/cstr.h"
#include "core/cstr_table.h"
#include "core/ensure.h"
#include "core/hashtable.h"
#include "core/ma.h"

struct CstrTable {
  Hashtable *strings;
};

static void free_cstr_table_entry(void *cstr_entry)
{
  ma_free(*(char**) cstr_entry);
}

CstrTable* cstr_table_new()
{
  HashElemInfo cstr_table = {
    ht_cstr_elem_hash, { free_cstr_table_entry }, sizeof (char*),
    ht_cstr_elem_cmp, NULL, NULL };
  CstrTable *table = ma_malloc(sizeof *table);
  table->strings = hashtable_new(cstr_table);
  return table;
}

void cstr_table_delete(CstrTable *table)
{
  if (!table) return;
  hashtable_delete(table->strings);
  ma_free(table);
}

void cstr_table_add(CstrTable *table, const char *cstr)
{
  char *dup;
  int rval;
  assert(table && cstr);
  assert(!cstr_table_get(table, cstr));
  dup = cstr_dup(cstr);
  rval = hashtable_add(table->strings, &dup);
  assert(rval == 1);
}

const char* cstr_table_get(const CstrTable *table, const char *cstr)
{
  const char **entry;
  assert(table && cstr);
  entry = hashtable_get(table->strings, &cstr);
  return entry ? *entry : NULL;
}

int cstr_table_unit_test(GT_Error *err)
{
  CstrTable *table;
  int had_err = 0;
  gt_error_check(err);
  table = cstr_table_new();
  ensure(had_err, !cstr_table_get(table, "foo"));
  if (!had_err)
    cstr_table_add(table, "foo");
  ensure(had_err, !strcmp(cstr_table_get(table, "foo"), "foo"));
  cstr_table_delete(table);
  return had_err;
}
