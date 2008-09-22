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

struct GtCstrTable {
  GtHashtable *strings;
};

static void free_gt_cstr_table_entry(void *cstr_entry)
{
  gt_free(*(char**) cstr_entry);
}

GtCstrTable* gt_cstr_table_new()
{
  HashElemInfo cstr_table = {
    gt_ht_cstr_elem_hash, { free_gt_cstr_table_entry }, sizeof (char*),
    gt_ht_cstr_elem_cmp, NULL, NULL };
  GtCstrTable *table = gt_malloc(sizeof *table);
  table->strings = gt_hashtable_new(cstr_table);
  return table;
}

void gt_cstr_table_delete(GtCstrTable *table)
{
  if (!table) return;
  gt_hashtable_delete(table->strings);
  gt_free(table);
}

void gt_cstr_table_add(GtCstrTable *table, const char *cstr)
{
  char *dup;
  int rval;
  assert(table && cstr);
  assert(!gt_cstr_table_get(table, cstr));
  dup = gt_cstr_dup(cstr);
  rval = gt_hashtable_add(table->strings, &dup);
  assert(rval == 1);
}

const char* gt_cstr_table_get(const GtCstrTable *table, const char *cstr)
{
  const char **entry;
  assert(table && cstr);
  entry = gt_hashtable_get(table->strings, &cstr);
  return entry ? *entry : NULL;
}

static enum iterator_op store_type(void *elem, void *data,
                                   GT_UNUSED GtError *err)
{
  GtStrArray *types = data;
  gt_error_check(err);
  assert(elem && types);
  gt_strarray_add_cstr(types, elem);
  return CONTINUE_ITERATION;
}

GtStrArray* gt_cstr_table_get_all(const GtCstrTable *table)
{
  int had_err;
  GtStrArray *cstrs;
  assert(table);
  cstrs = gt_strarray_new();
  had_err = gt_hashtable_foreach_ordered(table->strings, store_type, cstrs,
                                      (GtCompare) strcmp, NULL);
  assert(!had_err);
  return cstrs;
}

int gt_cstr_table_unit_test(GtError *err)
{
  GtCstrTable *table;
  int had_err = 0;
  gt_error_check(err);
  table = gt_cstr_table_new();
  ensure(had_err, !gt_cstr_table_get(table, "foo"));
  if (!had_err)
    gt_cstr_table_add(table, "foo");
  ensure(had_err, !strcmp(gt_cstr_table_get(table, "foo"), "foo"));
  gt_cstr_table_delete(table);
  return had_err;
}
