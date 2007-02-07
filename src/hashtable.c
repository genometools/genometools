/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
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
  info->iterfunc(key, value, info->data);
  return ST_CONTINUE;
}

void hashtable_foreach(Hashtable *ht, Hashiteratorfunc iterfunc, void *data)
{
  St_iterfunc_info info;
  assert(ht && iterfunc);
  info.iterfunc = iterfunc;
  info.data = data;
  (void) st_foreach(ht->st_table, st_iterfunc, (st_data_t) &info);
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

static void hashtable_test(Hash_type hash_type)
{
  char *s1 = "foo", *s2 = "bar";
  Hashtable *ht;

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
  assert(hashtable_get(ht, s1) == s2);
  assert(!hashtable_get(ht, s2));
  hashtable_free(ht);

  /* hashes containing two elements */
  ht = hashtable_new(hash_type, NULL, NULL);
  hashtable_add(ht, s1, s2);
  hashtable_add(ht, s2, s1);
  assert(hashtable_get(ht, s1) == s2);
  assert(hashtable_get(ht, s2) == s1);
  hashtable_free(ht);
}

int hashtable_unit_test(void)
{
  /* direct hash */
  hashtable_test(HASH_DIRECT);

  /* string hash */
  hashtable_test(HASH_STRING);

  return EXIT_SUCCESS;
}

void hashtable_free(Hashtable *ht)
{
  if (!ht) return;
  hashtable_reset(ht);
  assert(ht->st_table);
  st_free_table(ht->st_table);
  free(ht);
}
