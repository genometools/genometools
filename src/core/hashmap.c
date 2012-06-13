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

#include "core/cstr_api.h"
#include "core/hashtable.h"
#include "core/ma.h"
#include "core/hashmap.h"
#include "core/hashmap-generic.h"

/* Hashmaps are implemented as Hashtables */

struct map_entry
{
  void *key;
  void *value;
};

struct hm_freefuncs
{
  GtFree keyfree, valuefree;
};

static void
hm_elem_free(void *elem, void *table_data)
{
  struct hm_freefuncs *ff = table_data;
  gt_assert(elem && ff);
  if (ff->keyfree)
    ff->keyfree(((struct map_entry *)elem)->key);
  if (ff->valuefree)
    ff->valuefree(((struct map_entry *)elem)->value);
}

GtHashmap* gt_hashmap_new(GtHashType keyhashtype, GtFree keyfree,
                          GtFree valuefree)
{
  struct hm_freefuncs *ff = gt_malloc(sizeof (*ff));
  ff->keyfree = keyfree;
  ff->valuefree = valuefree;
  switch (keyhashtype) {
    case GT_HASH_DIRECT:
      {
        HashElemInfo hm_directkey_eleminfo = {
          gt_ht_ptr_elem_hash, { .free_elem_with_data = hm_elem_free },
          sizeof (struct map_entry), gt_ht_ptr_elem_cmp, ff, gt_free_func
        };
        return (GtHashmap*) gt_hashtable_new(hm_directkey_eleminfo);
      }
    case GT_HASH_STRING:
      {
        HashElemInfo hm_strkey_eleminfo = {
          gt_ht_cstr_elem_hash, { .free_elem_with_data = hm_elem_free },
          sizeof (struct map_entry), gt_ht_cstr_elem_cmp, ff, gt_free_func
        };
        return (GtHashmap*) gt_hashtable_new(hm_strkey_eleminfo);
      }
    default: gt_assert(0);
  }
  return NULL; /* for stupid compilers nagging too much */
}

GtHashmap* gt_hashmap_ref(GtHashmap *hm)
{
  return (GtHashmap*) gt_hashtable_ref((GtHashtable*) hm);
}

void* gt_hashmap_get(GtHashmap *hm, const void *key)
{
  struct map_entry *elem = gt_hashtable_get((GtHashtable*) hm, &key);
  return (elem!=NULL)?elem->value:NULL;
}

void* gt_hashmap_get_key(GtHashmap *hm, const void *key)
{
  struct map_entry *elem = gt_hashtable_get((GtHashtable*) hm, &key);
  return (elem!=NULL)?elem->key:NULL;
}

void gt_hashmap_add(GtHashmap *hm, void *key, void *value)
{
  struct map_entry keyvalpair = { key, value };
  if (!gt_hashtable_add((GtHashtable*) hm, &keyvalpair)) {
    ((struct map_entry *)gt_hashtable_get((GtHashtable*) hm, &keyvalpair))
      ->value = value;
  }
}

void gt_hashmap_remove(GtHashmap *hm, const void *key)
{
  gt_hashtable_remove((GtHashtable*) hm, &key);
}

/* iteration support structures and functions */
struct hashiteration_state
{
  GtHashmapVisitFunc visit;
  void *data;
  GtCompare keycmp;
};

static int
gt_hashmap_cmp(const void *elemA, const void *elemB, void *data)
{
  const struct map_entry *entryA = elemA,
    *entryB = elemB;
  struct hashiteration_state *dip = data;
  return ((GtCompareWithData)dip->keycmp)(entryA->key, entryB->key,
                                        dip->data);
}

static enum iterator_op
gt_hashmap_visit(void *elem, void *data, GtError *err)
{
  struct map_entry *me = elem;
  struct hashiteration_state *state = data;
  return state->visit(me->key, me->value, state->data, err) ?
    STOP_ITERATION : CONTINUE_ITERATION;
}

/* iterate over the hashmap in key order given by compare function <cmp> */
int gt_hashmap_foreach_ordered(GtHashmap *hm, GtHashmapVisitFunc visit,
                               void *data, GtCompare cmp, GtError *err)
{
  struct hashiteration_state state = { visit, data, cmp};
  return gt_hashtable_foreach_ordered((GtHashtable*) hm, gt_hashmap_visit,
                                      &state, (GtCompare) gt_hashmap_cmp, err);
}

int gt_hashmap_foreach(GtHashmap *hm, GtHashmapVisitFunc visit, void *data,
                       GtError *err)
{
  struct hashiteration_state state = { visit, data, NULL };
  return gt_hashtable_foreach((GtHashtable*) hm, gt_hashmap_visit, &state, err);
}

int gt_hashmap_foreach_in_key_order(GtHashmap *hm, GtHashmapVisitFunc iter,
                                    void *data, GtError *err)
{
  struct hashiteration_state state = { iter, data, NULL };
  return gt_hashtable_foreach_in_default_order((GtHashtable*) hm,
                                               gt_hashmap_visit, &state, err);
}

void gt_hashmap_reset(GtHashmap *hm)
{
  gt_hashtable_reset((GtHashtable*) hm);
}

#define my_ensure(err_state, predicate)         \
  if (!(predicate)) {                           \
    err_state = -1;                             \
    break;                                      \
  }

DECLARE_HASHMAP(unsigned long, testul, unsigned long long, testull,
                static, inline)
DEFINE_HASHMAP(unsigned long, testul, unsigned long long, testull,
               gt_ht_ul_elem_hash, gt_ht_ul_elem_cmp,
               NULL_DESTRUCTOR, NULL_DESTRUCTOR, static, inline)

static int
gt_hashmap_test(GtHashType hash_type)
{
  char *s1 = "foo", *s2 = "bar";
  GT_UNUSED unsigned long ul1 = 1UL, ul2 = 2UL;
  GT_UNUSED unsigned long long ull1 = 3ULL, ull2 = 4ULL, *sptr = NULL,
                               *tptr = NULL;
  GtHashmap *hm;
  GtHashtable *ht;
  int had_err = 0;
  do {
    /* empty hash */
    hm = gt_hashmap_new(hash_type, NULL, NULL);
    gt_hashmap_delete(hm);

    /* empty hash with reset */
    hm = gt_hashmap_new(hash_type, NULL, NULL);
    gt_hashmap_reset(hm);
    gt_hashmap_delete(hm);

    /* hashes containing one element */
    hm = gt_hashmap_new(hash_type, NULL, NULL);
    gt_hashmap_add(hm, s1, s2);
    my_ensure(had_err, gt_hashmap_get(hm, s1) == s2);
    my_ensure(had_err, !gt_hashmap_get(hm, s2));
    gt_hashmap_delete(hm);

    /* hashes containing two elements */
    hm = gt_hashmap_new(hash_type, NULL, NULL);
    gt_hashmap_add(hm, s1, s2);
    gt_hashmap_add(hm, s2, s1);
    my_ensure(had_err, gt_hashmap_get(hm, s1) == s2);
    my_ensure(had_err, gt_hashmap_get(hm, s2) == s1);

    /* remove element A and ensure it's no longer present */
    gt_hashmap_remove(hm, s1);
    my_ensure(had_err, !gt_hashmap_get(hm, s1));
    my_ensure(had_err, gt_hashmap_get(hm, s2) == s1);
    gt_hashmap_delete(hm);

    /* hashes containing two elements (store key and value in
     * hashmap) where simple free is the correct way to get rid of them
     */
    if (hash_type == GT_HASH_STRING)
    {
      hm = gt_hashmap_new(hash_type, gt_free_func, gt_free_func);

      gt_hashmap_add(hm, gt_cstr_dup(s1), gt_cstr_dup(s2));
      gt_hashmap_add(hm, gt_cstr_dup(s2), gt_cstr_dup(s1));
      my_ensure(had_err, !strcmp(gt_hashmap_get(hm, s1), s2));
      my_ensure(had_err, !strcmp(gt_hashmap_get(hm, s2), s1));
      gt_hashmap_remove(hm, s1); /* remove first element */
      my_ensure(had_err, !gt_hashmap_get(hm, s1));
      my_ensure(had_err, !strcmp(gt_hashmap_get(hm, s2),  s1));
      gt_hashmap_delete(hm);
    }

    /* test direct modification of hash contents */
    ht = testul_testull_gt_hashmap_new();
    my_ensure(had_err, !sptr);
    sptr = testul_testull_gt_hashmap_add_and_return_storage(ht, ul1, ull1);
    my_ensure(had_err, *sptr == ull1);
    sptr = testul_testull_gt_hashmap_add_and_return_storage(ht, ul1, ull2);
    my_ensure(had_err, *sptr == ull2);
    (*sptr)++;
    tptr = testul_testull_gt_hashmap_get(ht, ul1);
    my_ensure(had_err, tptr == sptr);
    my_ensure(had_err, *tptr == ull2+1);
    gt_hashtable_delete(ht);

  } while (0);
  return had_err;
}

int
gt_hashmap_unit_test(GtError *err)
{
  int had_err;
  gt_error_check(err);

  /* direct hash */
  had_err = gt_hashmap_test(GT_HASH_DIRECT);

  /* string hash */
  if (!had_err)
    had_err = gt_hashmap_test(GT_HASH_STRING);

  if (had_err)
  {
    gt_error_set(err, "hashmap operation created inconsistent state.");
  }
  return had_err;
}

void gt_hashmap_delete(GtHashmap *hm)
{
  gt_hashtable_delete((GtHashtable*) hm);
}
