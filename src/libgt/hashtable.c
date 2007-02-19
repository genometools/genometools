/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "ensure.h"
#include "hashtable.h"
#include "st.h"
#include "xansi.h"

struct Hashtable
{
  Hash_type hash_type;
  Hashkeyfreefunc key_free;
  Hashvaluefreefunc value_free;
  st_table *st_table;
};

typedef struct {
  Hashiteratorfunc iterfunc;
  void *data;
  Error *err;
  int has_err;
} St_iterfunc_info;

Hashtable* hashtable_new(Hash_type hash_type, Hashkeyfreefunc key_free,
                         Hashvaluefreefunc value_free)
{
  Hashtable *ht = xmalloc(sizeof(Hashtable));
  ht->hash_type = hash_type;
  ht->key_free = key_free;
  ht->value_free = value_free;
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
  assert(ht && key);
  (void) st_lookup(ht->st_table, (st_data_t) key, &value);
  vptr = (void*) value;
  return vptr;
}

void hashtable_add(Hashtable *ht, void *key, void *value)
{
  assert(ht && key && value);
  /* the hash table containes no element with this key yet */
  assert(!st_lookup(ht->st_table, (st_data_t) key, NULL));
  (void) st_insert(ht->st_table, (st_data_t) key, (st_data_t) value);
}

void hashtable_remove(Hashtable *ht, void *key)
{
  st_data_t value = 0;
  st_data_t key_t = (st_data_t) key;
  assert(ht && key);
  (void) st_lookup(ht->st_table, (st_data_t) key, &value);
  /* the hashtable containes an element with this key already */
  assert(value);
  (void) st_delete(ht->st_table, &key_t, NULL);
  if (ht->key_free)
    ht->key_free(key);
  if (ht->value_free)
    ht->value_free(&value);
}

static int st_iterfunc(void *key, void *value, void *data)
{
  St_iterfunc_info *info = (St_iterfunc_info*) data;
  assert(info->iterfunc);
  info->has_err = info->iterfunc(key, value, info->data, info->err);
  if (info->has_err)
    return ST_STOP;
  return ST_CONTINUE;
}

int hashtable_foreach(Hashtable *ht, Hashiteratorfunc iterfunc, void *data,
                      Error *err)
{
  St_iterfunc_info info;
  error_check(err);
  assert(ht && iterfunc);
  info.iterfunc = iterfunc;
  info.data = data;
  info.err = err;
  info.has_err = 0;
  (void) st_foreach(ht->st_table, st_iterfunc, (st_data_t) &info);
  return info.has_err;
}

static int remove_key_value_pair(void *key, void *value, void *data)
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
  (void) st_foreach(ht->st_table, remove_key_value_pair, (st_data_t) ht);
}

static int hashtable_test(Hash_type hash_type, Error *err)
{
  char *s1 = "foo", *s2 = "bar";
  Hashtable *ht;
  int has_err = 0;
  error_check(err);

  /* empty hash */
  ht = hashtable_new(hash_type, NULL, NULL);
  hashtable_free(ht);

  /* empty hash with reset */
  ht = hashtable_new(hash_type, NULL, NULL);
  hashtable_reset(ht);
  hashtable_free(ht);

  /* hashes containing one element */
  ht = hashtable_new(hash_type, NULL, NULL);
  hashtable_add(ht, s1, s2);
  ensure(has_err, hashtable_get(ht, s1) == s2);
  ensure(has_err, !hashtable_get(ht, s2));
  hashtable_free(ht);

  /* hashes containing two elements */
  ht = hashtable_new(hash_type, NULL, NULL);
  hashtable_add(ht, s1, s2);
  hashtable_add(ht, s2, s1);
  ensure(has_err, hashtable_get(ht, s1) == s2);
  ensure(has_err, hashtable_get(ht, s2) == s1);
  hashtable_free(ht);

  return has_err;
}

int hashtable_unit_test(Error *err)
{
  int has_err;
  error_check(err);

  /* direct hash */
  has_err = hashtable_test(HASH_DIRECT, err);

  /* string hash */
  if (!has_err)
    has_err = hashtable_test(HASH_STRING, err);

  return has_err;
}

void hashtable_free(Hashtable *ht)
{
  if (!ht) return;
  hashtable_reset(ht);
  assert(ht->st_table);
  st_free_table(ht->st_table);
  free(ht);
}
