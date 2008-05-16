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

#include <assert.h>
#include "libgtcore/array.h"
#include "libgtcore/cstr.h"
#include "libgtcore/ensure.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/st.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"

struct Hashtable
{
  HashType hash_type;
  FreeFunc key_free;
  FreeFunc value_free;
  st_table *st_table;
};

typedef struct {
  Hashiteratorfunc iterfunc;
  void *data;
  Error *err;
  int had_err;
} St_iterfunc_info;

Hashtable* hashtable_new(HashType hash_type, FreeFunc keyfree,
                         FreeFunc valuefree)
{
  Hashtable *ht = ma_malloc(sizeof (Hashtable));
  ht->hash_type = hash_type;
  ht->key_free = keyfree;
  ht->value_free = valuefree;
  switch (hash_type) {
    case HASH_DIRECT:
      ht->st_table = st_init_numtable();
      break;
    case HASH_STRING:
      ht->st_table = st_init_strtable();
  }
  assert(ht->st_table);
  return ht;
}

void* hashtable_get(Hashtable *ht, const void *key)
{
  st_data_t value = 0;
  void *vptr;
  assert(ht);
  (void) st_lookup(ht->st_table, (st_data_t) key, &value);
  vptr = (void*) value;
  return vptr;
}

void hashtable_add(Hashtable *ht, void *key, void *value)
{
  assert(ht && value);
  /* the hash table containes no element with this key yet */
  assert(!st_lookup(ht->st_table, (st_data_t) key, NULL));
  (void) st_insert(ht->st_table, (st_data_t) key, (st_data_t) value);
}

void hashtable_remove(Hashtable *ht, const void *key)
{
  st_data_t value = 0;
  st_data_t key_t = (st_data_t) key;
  assert(ht && key);
  (void) st_lookup(ht->st_table, (st_data_t) key, &value);
  /* the hashtable containes an element with this key already */
  assert(value);
  (void) st_delete(ht->st_table, &key_t, NULL);
  if (ht->key_free)
    ht->key_free((void*) key_t);
  if (ht->value_free)
    ht->value_free((void*) value);
}

static int st_iterfunc(void *key, void *value, void *data, UNUSED Error *err)
{
  St_iterfunc_info *info = (St_iterfunc_info*) data;
  assert(info->iterfunc);
  info->had_err = info->iterfunc(key, value, info->data, info->err);
  if (info->had_err)
    return ST_STOP;
  return ST_CONTINUE;
}

int hashtable_foreach(Hashtable *ht, Hashiteratorfunc iterfunc, void *data,
                      Error *err)
{
  St_iterfunc_info info;
  assert(ht && iterfunc);
  info.iterfunc = iterfunc;
  info.data = data;
  info.err = err;
  info.had_err = 0;
  (void) st_foreach(ht->st_table, st_iterfunc, (st_data_t) &info, err);
  return info.had_err;
}

typedef struct {
  void *key,
       *value;
} HashEntry;

static int save_hash_entry(void *key, void *value, void *data,
                           UNUSED Error *err)
{
  Array *hash_entries;
  HashEntry he;
  error_check(err);
  assert(value && data);
  hash_entries = (Array*) data;
  he.key = key;
  he.value = value;
  array_add(hash_entries, he);
  return 0;
}

static int ulongcmp(const void *a, const void *b)
{
  if ((unsigned long) a < (unsigned long) b)
    return -1;
  if ((unsigned long) a == (unsigned long) b)
    return 0;
  return 1;
}

/* XXX: remove this global variable, when qsort_r() has been written */
static Compare global_cmp = NULL;

static int compare_hash_entries(const void *a, const void *b)
{
  HashEntry *he_a = (HashEntry*) a, *he_b = (HashEntry*) b;
  assert(he_a && he_b && global_cmp);
  return global_cmp(he_a->key, he_b->key);
}

int hashtable_foreach_ordered(Hashtable *ht, Hashiteratorfunc iterfunc,
                              void *data, Compare cmp, Error *err)
{
  Array *hash_entries;
  HashEntry *he;
  unsigned long i;
  int had_err;
  assert(ht && iterfunc && cmp);
  hash_entries = array_new(sizeof (HashEntry));
  had_err = hashtable_foreach(ht, save_hash_entry, hash_entries, err);
  if (!had_err) {
    global_cmp = cmp;
    qsort(array_get_space(hash_entries), array_size(hash_entries),
          array_elem_size(hash_entries), compare_hash_entries);
    global_cmp = NULL;
    for (i = 0; !had_err && i < array_size(hash_entries); i++) {
      he = array_get(hash_entries, i);
      had_err = iterfunc(he->key, he->value, data, err);
    }
  }
  array_delete(hash_entries);
  return had_err;
}

int hashtable_foreach_ao(Hashtable *ht, Hashiteratorfunc iterfunc, void *data,
                         Error *err)
{
  assert(ht && iterfunc);
  assert(ht->hash_type == HASH_STRING);
  return hashtable_foreach_ordered(ht, iterfunc, data, (Compare) strcmp, err);
}

int hashtable_foreach_no(Hashtable *ht, Hashiteratorfunc iterfunc, void *data,
                         Error *err)
{
  assert(ht && iterfunc);
  assert(ht->hash_type == HASH_DIRECT);
  return hashtable_foreach_ordered(ht, iterfunc, data, ulongcmp, err);
}

static int remove_key_value_pair(void *key, void *value, void *data,
                                 UNUSED Error *err)
{
  Hashtable *ht= (Hashtable*) data;
  assert(ht);
  if (ht->key_free)
    ht->key_free(key);
  if (ht->value_free)
    ht->value_free(value);
  return ST_DELETE;
}

void hashtable_reset(Hashtable *ht)
{
  assert(ht && ht->st_table);
  (void) st_foreach(ht->st_table, remove_key_value_pair, (st_data_t) ht, NULL);
}

static int hashtable_test(HashType hash_type, Error *err)
{
  char *s1 = "foo", *s2 = "bar";
  Hashtable *ht;
  int had_err = 0;
  error_check(err);

  /* empty hash */
  ht = hashtable_new(hash_type, NULL, NULL);
  hashtable_delete(ht);

  /* empty hash with reset */
  ht = hashtable_new(hash_type, NULL, NULL);
  hashtable_reset(ht);
  hashtable_delete(ht);

  /* hashes containing one element */
  ht = hashtable_new(hash_type, NULL, NULL);
  hashtable_add(ht, s1, s2);
  ensure(had_err, hashtable_get(ht, s1) == s2);
  ensure(had_err, !hashtable_get(ht, s2));
  hashtable_delete(ht);

  /* hashes containing two elements */
  if (!had_err) {
    ht = hashtable_new(hash_type, NULL, NULL);
    hashtable_add(ht, s1, s2);
    hashtable_add(ht, s2, s1);
    ensure(had_err, hashtable_get(ht, s1) == s2);
    ensure(had_err, hashtable_get(ht, s2) == s1);
    hashtable_remove(ht, s1); /* remove first element */
    ensure(had_err, !hashtable_get(ht, s1));
    ensure(had_err, hashtable_get(ht, s2) == s1);
    hashtable_delete(ht);
  }

  /* hashes containing two elements (store key and value in hashtable) */
  if (!had_err && hash_type == HASH_STRING) {
    ht = hashtable_new(hash_type, ma_free_func, ma_free_func);
    hashtable_add(ht, cstr_dup(s1), cstr_dup(s2));
    hashtable_add(ht, cstr_dup(s2), cstr_dup(s1));
    ensure(had_err, !strcmp(hashtable_get(ht, s1), s2));
    ensure(had_err, !strcmp(hashtable_get(ht, s2), s1));
    hashtable_remove(ht, s1); /* remove first element */
    ensure(had_err, !hashtable_get(ht, s1));
    ensure(had_err, !strcmp(hashtable_get(ht, s2),  s1));
    hashtable_delete(ht);
  }

  return had_err;
}

int hashtable_unit_test(Error *err)
{
  int had_err;
  error_check(err);

  /* direct hash */
  had_err = hashtable_test(HASH_DIRECT, err);

  /* string hash */
  if (!had_err)
    had_err = hashtable_test(HASH_STRING, err);

  return had_err;
}

void hashtable_delete(Hashtable *ht)
{
  if (!ht) return;
  hashtable_reset(ht);
  assert(ht->st_table);
  st_free_table(ht->st_table);
  ma_free(ht);
}
