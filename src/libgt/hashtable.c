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
  FreeFunc key_free;
  FreeFunc value_free;
  st_table *st_table;
};

typedef struct {
  Hashiteratorfunc iterfunc;
  void *data;
  Env *env;
  int has_err;
} St_iterfunc_info;

Hashtable* hashtable_new(Hash_type hash_type, FreeFunc keyfree,
                         FreeFunc valuefree, Env *env)
{
  Hashtable *ht = env_ma_malloc(env, sizeof (Hashtable));
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

void hashtable_remove(Hashtable *ht, void *key, Env *env)
{
  st_data_t value = 0;
  st_data_t key_t = (st_data_t) key;
  assert(ht && key);
  (void) st_lookup(ht->st_table, (st_data_t) key, &value);
  /* the hashtable containes an element with this key already */
  assert(value);
  (void) st_delete(ht->st_table, &key_t, NULL);
  if (ht->key_free)
    ht->key_free(key, env);
  if (ht->value_free)
    ht->value_free(&value, env);
}

static int st_iterfunc(void *key, void *value, void *data)
{
  St_iterfunc_info *info = (St_iterfunc_info*) data;
  assert(info->iterfunc);
  info->has_err = info->iterfunc(key, value, info->data, info->env);
  if (info->has_err)
    return ST_STOP;
  return ST_CONTINUE;
}

int hashtable_foreach(Hashtable *ht, Hashiteratorfunc iterfunc, void *data,
                      Env *env)
{
  St_iterfunc_info info;
  env_error_check(env);
  assert(ht && iterfunc);
  info.iterfunc = iterfunc;
  info.data = data;
  info.env = env;
  info.has_err = 0;
  (void) st_foreach(ht->st_table, st_iterfunc, (st_data_t) &info);
  return info.has_err;
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

void hashtable_reset(Hashtable *ht)
{
  assert(ht && ht->st_table);
  (void) st_foreach(ht->st_table, remove_key_value_pair, (st_data_t) ht);
}

static int hashtable_test(Hash_type hash_type, Env *env)
{
  char *s1 = "foo", *s2 = "bar";
  Hashtable *ht;
  int has_err = 0;
  env_error_check(env);

  /* empty hash */
  ht = hashtable_new(hash_type, NULL, NULL, env);
  hashtable_delete(ht, env);

  /* empty hash with reset */
  ht = hashtable_new(hash_type, NULL, NULL, env);
  hashtable_reset(ht);
  hashtable_delete(ht, env);

  /* hashes containing one element */
  ht = hashtable_new(hash_type, NULL, NULL, env);
  hashtable_add(ht, s1, s2);
  ensure(has_err, hashtable_get(ht, s1) == s2);
  ensure(has_err, !hashtable_get(ht, s2));
  hashtable_delete(ht, env);

  /* hashes containing two elements */
  ht = hashtable_new(hash_type, NULL, NULL, env);
  hashtable_add(ht, s1, s2);
  hashtable_add(ht, s2, s1);
  ensure(has_err, hashtable_get(ht, s1) == s2);
  ensure(has_err, hashtable_get(ht, s2) == s1);
  hashtable_delete(ht, env);

  return has_err;
}

int hashtable_unit_test(Env *env)
{
  int has_err;
  env_error_check(env);

  /* direct hash */
  has_err = hashtable_test(HASH_DIRECT, env);

  /* string hash */
  if (!has_err)
    has_err = hashtable_test(HASH_STRING, env);

  return has_err;
}

void hashtable_delete(Hashtable *ht, Env *env)
{
  if (!ht) return;
  hashtable_reset(ht);
  assert(ht->st_table);
  st_free_table(ht->st_table);
  env_ma_free(ht, env);
}
