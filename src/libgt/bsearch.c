/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "bsearch.h"
#include "ensure.h"

void bsearch_all(Array *members, const void *key, const void *base,
                 size_t nmemb, size_t size,
                 int (*compar) (const void *, const void *))
{
  bsearch_all_mark(members, key, base, nmemb, size, compar, NULL);
}

void bsearch_all_mark(Array *members, const void *key, const void *base,
                      size_t nmemb, size_t size,
                      int (*compar) (const void *, const void *), Bittab *b)
{
  const void *baseptr = base;
  const void *ptr, /* the current element we consider */
             *tmp_ptr;
  int limit, rval;

  assert(members && key && size && compar);
  assert(!b || bittab_size(b) == nmemb);

  /* the actual binary search */
  for (limit = nmemb; limit != 0; limit >>= 1) {
    ptr = baseptr + (limit >> 1) * size;
    if ((rval = compar(key, ptr)) == 0) {
      /* element found */
      array_add(members, ptr);
      if (b)
        bittab_set_bit(b, (ptr - base) / size); /* mark found element */
      /* looking left for equal elements */
      for (tmp_ptr = ptr - size;
           tmp_ptr >= base && !compar(key, tmp_ptr);
           tmp_ptr -= size) {
        array_add(members, tmp_ptr);
        if (b)
          bittab_set_bit(b, (tmp_ptr - base) / size); /* mark found element */
      }
      /* looking right for equal elements */
      for (tmp_ptr = ptr + size;
           tmp_ptr < base + nmemb * size && !compar(key, tmp_ptr);
           tmp_ptr += size) {
        array_add(members, tmp_ptr);
        if (b)
          bittab_set_bit(b, (tmp_ptr - base) / size); /* mark found element */
      }
      return;
    }
    if (rval > 0) {
      baseptr = ptr + size;
      limit--;
    }
  }
}

static int cmp(const void *a_ptr, const void *b_ptr)
{
  int a, b;
  assert(a_ptr && b_ptr);
  a = *(int*) a_ptr;
  b = *(int*) b_ptr;
  if (a == b)
    return 0;
  if (a < b)
    return -1;
  return 1;
}

/* XXX: This unit test could be done much better by filling an array randomly,
   sorting it, and comparing bsearch_all() against a brute force implementation.
*/
int bsearch_unit_test(Error *err)
{
  Array *elements, *members;
  int key, element, *member_ptr;
  Bittab *b;
  int has_err = 0;

  error_check(err);

  elements = array_new(sizeof(int));
  members = array_new(sizeof(int*));

  /* the empty case */
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof(int), cmp);
  ensure(has_err, !array_size(members)); /* no member found */

  /* 1 element */
  key = 7;
  element = 7;
  array_add(elements, element);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof(int), cmp);
  ensure(has_err, array_size(members) == 1); /* one member found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(has_err, *member_ptr == element);

  key = -7;
  array_set_size(members, 0);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof(int), cmp);
  ensure(has_err, !array_size(members)); /* no member found */

  /* 2 elements */
  key = 7;
  array_set_size(members, 0);
  array_add(elements, element);
  ensure(has_err, array_size(elements) == 2);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof(int), cmp);
  ensure(has_err, array_size(members) == 2); /* two members found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(has_err, *member_ptr == element);
  member_ptr = *(int**) array_get(members, 1);
  ensure(has_err, *member_ptr == element);

  key = -7;
  array_set_size(members, 0);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof(int), cmp);
  ensure(has_err, !array_size(members)); /* no member found */

  /* 3 elements */
  key = 7;
  array_set_size(members, 0);
  array_add(elements, element);
  ensure(has_err, array_size(elements) == 3);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof(int), cmp);
  ensure(has_err, array_size(members) == 3); /* three members found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(has_err, *member_ptr == element);
  member_ptr = *(int**) array_get(members, 1);
  ensure(has_err, *member_ptr == element);
  member_ptr = *(int**) array_get(members, 2);
  ensure(has_err, *member_ptr == element);

  key = -7;
  array_set_size(members, 0);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof(int), cmp);
  ensure(has_err, !array_size(members)); /* no member found */

  /* large case: -10 -5 -3 -3 -3 0 1 2 3 */
  array_set_size(elements, 0);
  element = -10;
  array_add(elements, element);
  element = -5;
  array_add(elements, element);
  element = -3;
  array_add(elements, element);
  array_add(elements, element);
  array_add(elements, element);
  element = 0;
  array_add(elements, element);
  element = 1;
  array_add(elements, element);
  element = 2;
  array_add(elements, element);
  element = 3;
  array_add(elements, element);
  ensure(has_err, array_size(elements) == 9);
  key = -3;
  array_set_size(members, 0);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof(int), cmp);
  ensure(has_err, array_size(members) == 3); /* three members found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(has_err, *member_ptr == -3);
  member_ptr = *(int**) array_get(members, 1);
  ensure(has_err, *member_ptr == -3);
  member_ptr = *(int**) array_get(members, 2);
  ensure(has_err, *member_ptr == -3);

  /* test bsearch_all_mark() with large case */
  array_set_size(members, 0);
  b = bittab_new(array_size(elements));
  bsearch_all_mark(members, &key, array_get_space(elements),
                   array_size(elements), sizeof(int), cmp, b);
  ensure(has_err, array_size(members) == 3); /* three members found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(has_err, *member_ptr == -3);
  member_ptr = *(int**) array_get(members, 1);
  ensure(has_err, *member_ptr == -3);
  member_ptr = *(int**) array_get(members, 2);
  ensure(has_err, *member_ptr == -3);
  /* the correct elements are marked (and only these) */
  ensure(has_err, bittab_bit_is_set(b, 2));
  ensure(has_err, bittab_bit_is_set(b, 3));
  ensure(has_err, bittab_bit_is_set(b, 4));
  ensure(has_err, bittab_count_set_bits(b) == 3);

  /* free */
  array_free(elements);
  array_free(members);
  bittab_free(b);

  return has_err;
}
