/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/ensure.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/st.h"
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
  Env *env;
  int had_err;
} St_iterfunc_info;

Hashtable* hashtable_new(HashType hash_type, FreeFunc keyfree,
                         FreeFunc valuefree, Env *env)
{
  Hashtable *ht = ma_malloc(sizeof (Hashtable));
  ht->hash_type = hash_type;
  ht->key_free = keyfree;
  ht->value_free = valuefree;
  switch (hash_type) {
    case HASH_DIRECT:
      ht->st_table = st_init_numtable(env);
      break;
    case HASH_STRING:
      ht->st_table = st_init_strtable(env);
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

void hashtable_add(Hashtable *ht, void *key, void *value, Env *env)
{
  assert(ht && value);
  /* the hash table containes no element with this key yet */
  assert(!st_lookup(ht->st_table, (st_data_t) key, NULL));
  (void) st_insert(ht->st_table, (st_data_t) key, (st_data_t) value, env);
}

void hashtable_remove(Hashtable *ht, void *key, Env *env)
{
  st_data_t value = 0;
  st_data_t key_t = (st_data_t) key;
  assert(ht && key);
  /*printf("hashtable_remove\n");*/
  (void) st_lookup(ht->st_table, (st_data_t) key, &value);
  /* the hashtable containes an element with this key already */
  assert(value);
  (void) st_delete(ht->st_table, &key_t, NULL, env);
  if (ht->key_free)
    ht->key_free(key, env);
  if (ht->value_free)
    ht->value_free((void*) value, env);
}

static int st_iterfunc(void *key, void *value, void *data,
                       /*@unused@*/ Env *env)
{
  St_iterfunc_info *info = (St_iterfunc_info*) data;
  assert(info->iterfunc);
  info->had_err = info->iterfunc(key, value, info->data, info->env);
  if (info->had_err)
    return ST_STOP;
  return ST_CONTINUE;
}

int hashtable_foreach(Hashtable *ht, Hashiteratorfunc iterfunc, void *data,
                      Env *env)
{
  St_iterfunc_info info;
  assert(ht && iterfunc);
  info.iterfunc = iterfunc;
  info.data = data;
  info.env = env;
  info.had_err = 0;
  (void) st_foreach(ht->st_table, st_iterfunc, (st_data_t) &info, env);
  return info.had_err;
}

typedef struct {
  void *key,
       *value;
} HashEntry;

static int save_hash_entry(void *key, void *value, void *data, Env *env)
{
  Array *hash_entries;
  HashEntry he;
  env_error_check(env);
  assert(key && value && data);
  hash_entries = (Array*) data;
  he.key = key;
  he.value = value;
  array_add(hash_entries, he, env);
  return 0;
}

int compare_hash_entries_alphabetically(const void *a, const void *b)
{
  HashEntry *he_a = (HashEntry*) a, *he_b = (HashEntry*) b;
  assert(he_a && he_b);
  return strcmp(he_a->key, he_b->key);
}

int compare_hash_entries_numerically(const void *a, const void *b)
{
  HashEntry *he_a = (HashEntry*) a, *he_b = (HashEntry*) b;
  assert(he_a && he_b);
  if ((unsigned long) he_a->key < (unsigned long) he_b->key)
    return -1;
  if ((unsigned long) he_a->key == (unsigned long) he_b->key)
    return 0;
  return 1;
}

int hashtable_foreach_ordered(Hashtable *ht, Hashiteratorfunc iterfunc,
                              void *data, int(*cmp)(const void*, const void*),
                              Env *env)
{
  Array *hash_entries;
  HashEntry *he;
  unsigned long i;
  int had_err;
  assert(ht && iterfunc && cmp);
  hash_entries = array_new(sizeof (HashEntry), env);
  had_err = hashtable_foreach(ht, save_hash_entry, hash_entries, env);
  if (!had_err) {
    qsort(array_get_space(hash_entries), array_size(hash_entries),
          array_elem_size(hash_entries), cmp);
    for (i = 0; !had_err && i < array_size(hash_entries); i++) {
      he = array_get(hash_entries, i);
      had_err = iterfunc(he->key, he->value, data, env);
    }
  }
  array_delete(hash_entries, env);
  return had_err;
}

int hashtable_foreach_ao(Hashtable *ht, Hashiteratorfunc iterfunc, void *data,
                         Env *env)
{
  assert(ht && iterfunc);
  assert(ht->hash_type == HASH_STRING);
  return hashtable_foreach_ordered(ht, iterfunc, data,
                                   compare_hash_entries_alphabetically, env);
}

int hashtable_foreach_no(Hashtable *ht, Hashiteratorfunc iterfunc, void *data,
                         Env *env)
{
  assert(ht && iterfunc);
  assert(ht->hash_type == HASH_DIRECT);
  return hashtable_foreach_ordered(ht, iterfunc, data,
                                   compare_hash_entries_numerically, env);
}

static int remove_key_value_pair(void *key, void *value, void *data, Env *env)
{
  Hashtable *ht= (Hashtable*) data;
  assert(ht);
  if (ht->key_free)
    ht->key_free(key, env);
  if (ht->value_free)
    ht->value_free(value, env);
  return ST_DELETE;
}

void hashtable_reset(Hashtable *ht, Env *env)
{
  assert(ht && ht->st_table);
  (void) st_foreach(ht->st_table, remove_key_value_pair, (st_data_t) ht, env);
}

static int hashtable_test(HashType hash_type, Env *env)
{
  char *s1 = "foo", *s2 = "bar";
  Hashtable *ht;
  int had_err = 0;
  env_error_check(env);

  /* empty hash */
  ht = hashtable_new(hash_type, NULL, NULL, env);
  hashtable_delete(ht, env);

  /* empty hash with reset */
  ht = hashtable_new(hash_type, NULL, NULL, env);
  hashtable_reset(ht, env);
  hashtable_delete(ht, env);

  /* hashes containing one element */
  ht = hashtable_new(hash_type, NULL, NULL, env);
  hashtable_add(ht, s1, s2, env);
  ensure(had_err, hashtable_get(ht, s1) == s2);
  ensure(had_err, !hashtable_get(ht, s2));
  hashtable_delete(ht, env);

  /* hashes containing two elements */
  ht = hashtable_new(hash_type, NULL, NULL, env);
  hashtable_add(ht, s1, s2, env);
  hashtable_add(ht, s2, s1, env);
  ensure(had_err, hashtable_get(ht, s1) == s2);
  ensure(had_err, hashtable_get(ht, s2) == s1);
  hashtable_delete(ht, env);

  return had_err;
}

int hashtable_unit_test(Env *env)
{
  int had_err;
  env_error_check(env);

  /* direct hash */
  had_err = hashtable_test(HASH_DIRECT, env);

  /* string hash */
  if (!had_err)
    had_err = hashtable_test(HASH_STRING, env);

  return had_err;
}

void hashtable_delete(Hashtable *ht, Env *env)
{
  if (!ht) return;
  hashtable_reset(ht, env);
  assert(ht->st_table);
  st_free_table(ht->st_table, env);
  ma_free(ht);
}
