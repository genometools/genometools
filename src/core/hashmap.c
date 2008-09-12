/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>
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
#include "core/hashtable.h"
#include "core/ma.h"
#define Hashmap Hashtable
#include "core/hashmap.h"

struct map_entry
{
  void *key;
  void *value;
};

struct hm_freefuncs
{
  GT_FreeFunc keyfree, valuefree;
};

static void
hm_elem_free(void *elem, void *table_data)
{
  struct hm_freefuncs *ff = table_data;
  assert(elem && ff);
  if (ff->keyfree)
    ff->keyfree(((struct map_entry *)elem)->key);
  if (ff->valuefree)
    ff->valuefree(((struct map_entry *)elem)->value);
}

extern Hashmap *
hashmap_new(HashType keyhashtype, GT_FreeFunc keyfree, GT_FreeFunc valuefree)
{
  struct hm_freefuncs *ff = gt_malloc(sizeof (*ff));
  ff->keyfree = keyfree;
  ff->valuefree = valuefree;
  switch (keyhashtype)
  {
  case HASH_DIRECT:
    {
      HashElemInfo hm_directkey_eleminfo = {
        ht_ptr_elem_hash, { .free_elem_with_data = hm_elem_free },
        sizeof (struct map_entry), ht_ptr_elem_cmp, ff, gt_free_func
      };
      return hashtable_new(hm_directkey_eleminfo);
    }
  case HASH_STRING:
    {
      HashElemInfo hm_strkey_eleminfo = {
        ht_cstr_elem_hash, { .free_elem_with_data = hm_elem_free },
        sizeof (struct map_entry), ht_cstr_elem_cmp, ff, gt_free_func
      };
      return hashtable_new(hm_strkey_eleminfo);
    }
  }
  fprintf(stderr, "Illegal hashtype requested!");
  abort();
  return NULL;                  /* for stupid compilers nagging too much */
}

extern void *
hashmap_get(Hashmap *hm, const void *key)
{
  struct map_entry *elem = hashtable_get(hm, &key);
  return (elem!=NULL)?elem->value:NULL;
}

extern void
hashmap_add(Hashmap *hm, void *key, void *value)
{
  struct map_entry keyvalpair = { key, value };
  if (hashtable_add(hm, &keyvalpair))
    ;
  else
    ((struct map_entry *)hashtable_get(hm, &keyvalpair))->value = value;
}

extern void
hashmap_remove(Hashmap *hm, const void *key)
{
  hashtable_remove(hm, &key);
}

/* iteration support structures and functions */
struct hashiteration_state
{
  Mapentryvisitfunc visit;
  void *data;
  GtCompare keycmp;
};

static int
hashmap_cmp(const void *elemA, const void *elemB, void *data)
{
  const struct map_entry *entryA = elemA,
    *entryB = elemB;
  struct hashiteration_state *dip = data;
  return ((GtCompareWithData)dip->keycmp)(entryA->key, entryB->key,
                                        dip->data);
}

static enum iterator_op
hashmap_visit(void *elem, void *data, GtError *err)
{
  struct map_entry *me = elem;
  struct hashiteration_state *state = data;
  return state->visit(me->key, me->value, state->data, err) ?
    STOP_ITERATION : CONTINUE_ITERATION;
}

/* iterate over the hashmap in key order given by compare function <cmp> */
extern int
hashmap_foreach_ordered(Hashmap *hm, Mapentryvisitfunc visit, void *data,
                        GtCompare cmp, GtError *err)
{
  struct hashiteration_state state = { visit, data, cmp};
  return hashtable_foreach_ordered(hm, hashmap_visit, &state,
                                   (GtCompare)hashmap_cmp, err);
}

extern int
hashmap_foreach(Hashmap *hm, Mapentryvisitfunc visit, void *data, GtError *err)
{
  struct hashiteration_state state = { visit, data, NULL };
  return hashtable_foreach(hm, hashmap_visit, &state, err);
}

extern int
hashmap_foreach_in_key_order(Hashmap *hm, Mapentryvisitfunc iter,
                             void *data, GtError *err)
{
  struct hashiteration_state state = { iter, data, NULL };
  return hashtable_foreach_in_default_order(hm, hashmap_visit, &state, err);
}

extern void
hashmap_reset(Hashmap *hm)
{
  hashtable_reset(hm);
}

#define my_ensure(err_state, predicate)         \
  if (!(predicate)) {                           \
    err_state = -1;                             \
    break;                                      \
  }

static int
hashmap_test(HashType hash_type)
{
  char *s1 = "foo", *s2 = "bar";
  Hashmap *hm;
  int had_err = 0;
  do {
    /* empty hash */
    hm = hashmap_new(hash_type, NULL, NULL);
    hashmap_delete(hm);

    /* empty hash with reset */
    hm = hashmap_new(hash_type, NULL, NULL);
    hashmap_reset(hm);
    hashmap_delete(hm);

    /* hashes containing one element */
    hm = hashmap_new(hash_type, NULL, NULL);
    hashmap_add(hm, s1, s2);
    my_ensure(had_err, hashmap_get(hm, s1) == s2);
    my_ensure(had_err, !hashmap_get(hm, s2));
    hashmap_delete(hm);

    /* hashes containing two elements */
    hm = hashmap_new(hash_type, NULL, NULL);
    hashmap_add(hm, s1, s2);
    hashmap_add(hm, s2, s1);
    my_ensure(had_err, hashmap_get(hm, s1) == s2);
    my_ensure(had_err, hashmap_get(hm, s2) == s1);

    /* remove element A and ensure it's no longer present */
    hashmap_remove(hm, s1);
    my_ensure(had_err, !hashmap_get(hm, s1));
    my_ensure(had_err, hashmap_get(hm, s2) == s1);
    hashmap_delete(hm);

    /* hashes containing two elements (store key and value in
     * hashmap) where simple free is the correct way to get rid of them
     */
    if (hash_type == HASH_STRING)
    {
      hm = hashmap_new(hash_type, gt_free_func, gt_free_func);

      hashmap_add(hm, gt_cstr_dup(s1), gt_cstr_dup(s2));
      hashmap_add(hm, gt_cstr_dup(s2), gt_cstr_dup(s1));
      my_ensure(had_err, !strcmp(hashmap_get(hm, s1), s2));
      my_ensure(had_err, !strcmp(hashmap_get(hm, s2), s1));
      hashmap_remove(hm, s1); /* remove first element */
      my_ensure(had_err, !hashmap_get(hm, s1));
      my_ensure(had_err, !strcmp(hashmap_get(hm, s2),  s1));
      hashmap_delete(hm);
    }
  } while (0);
  return had_err;
}

int
hashmap_unit_test(GtError *err)
{
  int had_err;
  gt_error_check(err);

  /* direct hash */
  had_err = hashmap_test(HASH_DIRECT);

  /* string hash */
  if (!had_err)
    had_err = hashmap_test(HASH_STRING);

  if (had_err)
  {
    gt_error_set(err, "hashmap operation created inconsistent state.");
  }
  return had_err;
}

extern void
hashmap_delete(Hashmap *hm)
{
  hashtable_delete(hm);
}
