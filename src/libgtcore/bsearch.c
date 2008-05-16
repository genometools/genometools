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
#include "libgtcore/bsearch.h"
#include "libgtcore/ensure.h"
#include "libgtcore/unused.h"

static void* bsearch_generic(Array *members, const void *key, const void *base,
                             size_t nmemb, size_t size, Compar compar,
                             void *data, Bittab *b)
{
  char *baseptr = (char *)base, *tmp_ptr,
   *ptr; /* the current element we consider */
  int limit, rval;

  assert(key && size && compar);
  assert(!b || bittab_size(b) == nmemb);

  /* the actual binary search */
  for (limit = nmemb; limit != 0; limit >>= 1) {
    ptr = baseptr + (limit >> 1) * size;
    if ((rval = compar(key, ptr, data)) == 0) {
      /* element found */
      if (members)
        array_add(members, ptr);
      else
        return ptr; /* because we got no <members> array */
      if (b)
        bittab_set_bit(b, (ptr - (char *)base) / size); /* mark found element */
      /* looking left for equal elements */
      for (tmp_ptr = ptr - size;
           tmp_ptr >= (char *)base && !compar(key, tmp_ptr, data);
           tmp_ptr -= size) {
        array_add(members, tmp_ptr);
        if (b)
          bittab_set_bit(
            b, (tmp_ptr - (char *)base) / size); /* mark found element */
      }
      /* looking right for equal elements */
      for (tmp_ptr = ptr + size;
           tmp_ptr <(char *) base + nmemb * size && !compar(key, tmp_ptr, data);
           tmp_ptr += size) {
        array_add(members, tmp_ptr);
        if (b)
          bittab_set_bit(
            b, (tmp_ptr - (char *)base) / size); /* mark found element */
      }
      return ptr;
    }
    if (rval > 0) {
      baseptr = ptr + size;
      limit--;
    }
  }
  return NULL;
}

void* bsearch_data(const void *key, const void *base, size_t nmemb, size_t size,
                   Compar compar, void *data)
{
  return bsearch_generic(NULL, key, base, nmemb, size, compar, data, NULL);
}

void bsearch_all(Array *members, const void *key, const void *base,
                 size_t nmemb, size_t size, Compar compar, void *data)
{
  assert(members);
  bsearch_generic(members, key, base, nmemb, size, compar, data, NULL);
}

void bsearch_all_mark(Array *members, const void *key, const void *base,
                      size_t nmemb, size_t size, Compar compar, void *data,
                      Bittab *b)
{
  assert(members);
  bsearch_generic(members, key, base, nmemb, size, compar, data, b);
}

static int cmp(const void *a_ptr, const void *b_ptr, UNUSED const void *unsused)
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
  int had_err = 0;

  error_check(err);

  elements = array_new(sizeof (int));
  members = array_new(sizeof (int*));

  /* the empty case */
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof (int), cmp, NULL);
  ensure(had_err, !array_size(members)); /* no member found */
  ensure(had_err, !bsearch_data(&key, array_get_space(elements),
                                array_size(elements), sizeof (int), cmp, NULL));

  /* 1 element */
  key = 7;
  element = 7;
  array_add(elements, element);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof (int), cmp, NULL);
  ensure(had_err, array_size(members) == 1); /* one member found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(had_err, *member_ptr == element);
  member_ptr = bsearch_data(&key, array_get_space(elements),
                            array_size(elements), sizeof (int), cmp, NULL);
  ensure(had_err, *member_ptr == element);

  key = -7;
  array_reset(members);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof (int), cmp, NULL);
  ensure(had_err, !array_size(members)); /* no member found */

  /* 2 elements */
  key = 7;
  array_reset(members);
  array_add(elements, element);
  ensure(had_err, array_size(elements) == 2);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof (int), cmp, NULL);
  ensure(had_err, array_size(members) == 2); /* two members found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(had_err, *member_ptr == element);
  member_ptr = *(int**) array_get(members, 1);
  ensure(had_err, *member_ptr == element);
  ensure(had_err, bsearch_data(&key, array_get_space(elements),
                               array_size(elements), sizeof (int), cmp, NULL));

  key = -7;
  array_reset(members);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof (int), cmp, NULL);
  ensure(had_err, !array_size(members)); /* no member found */
  ensure(had_err, !bsearch_data(&key, array_get_space(elements),
                                array_size(elements), sizeof (int), cmp, NULL));

  /* 3 elements */
  key = 7;
  array_reset(members);
  array_add(elements, element);
  ensure(had_err, array_size(elements) == 3);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof (int), cmp, NULL);
  ensure(had_err, array_size(members) == 3); /* three members found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(had_err, *member_ptr == element);
  member_ptr = *(int**) array_get(members, 1);
  ensure(had_err, *member_ptr == element);
  member_ptr = *(int**) array_get(members, 2);
  ensure(had_err, *member_ptr == element);
  ensure(had_err, bsearch_data(&key, array_get_space(elements),
                               array_size(elements), sizeof (int), cmp, NULL));

  key = -7;
  array_reset(members);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof (int), cmp, NULL);
  ensure(had_err, !array_size(members)); /* no member found */
  ensure(had_err, !bsearch_data(&key, array_get_space(elements),
                                array_size(elements), sizeof (int), cmp, NULL));

  /* large case: -10 -5 -3 -3 -3 0 1 2 3 */
  array_reset(elements);
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
  ensure(had_err, array_size(elements) == 9);
  key = -3;
  array_reset(members);
  bsearch_all(members, &key, array_get_space(elements), array_size(elements),
              sizeof (int), cmp, NULL);
  ensure(had_err, array_size(members) == 3); /* three members found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(had_err, *member_ptr == -3);
  member_ptr = *(int**) array_get(members, 1);
  ensure(had_err, *member_ptr == -3);
  member_ptr = *(int**) array_get(members, 2);
  ensure(had_err, *member_ptr == -3);
  ensure(had_err, bsearch_data(&key, array_get_space(elements),
                               array_size(elements), sizeof (int), cmp, NULL));

  /* test bsearch_all_mark() with large case */
  array_reset(members);
  b = bittab_new(array_size(elements));
  bsearch_all_mark(members, &key, array_get_space(elements),
                   array_size(elements), sizeof (int), cmp, NULL, b);
  ensure(had_err, array_size(members) == 3); /* three members found */
  member_ptr = *(int**) array_get(members, 0);
  ensure(had_err, *member_ptr == -3);
  member_ptr = *(int**) array_get(members, 1);
  ensure(had_err, *member_ptr == -3);
  member_ptr = *(int**) array_get(members, 2);
  ensure(had_err, *member_ptr == -3);
  /* the correct elements are marked (and only these) */
  ensure(had_err, bittab_bit_is_set(b, 2));
  ensure(had_err, bittab_bit_is_set(b, 3));
  ensure(had_err, bittab_bit_is_set(b, 4));
  ensure(had_err, bittab_count_set_bits(b) == 3);

  /* free */
  array_delete(elements);
  array_delete(members);
  bittab_delete(b);

  return had_err;
}
