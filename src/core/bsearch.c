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

#include "core/assert_api.h"
#include "core/bsearch.h"
#include "core/ensure.h"
#include "core/unused_api.h"

static void* bsearch_generic(GtArray *members, const void *key,
                             const void *base, size_t nmemb, size_t size,
                             GtCompareWithData compar, void *data,
                             GtBittab *b)
{
  char *baseptr = (char *)base, *tmp_ptr,
   *ptr; /* the current element we consider */
  int limit, rval;

  gt_assert(key && size && compar);
  gt_assert(!b || gt_bittab_size(b) == nmemb);

  /* the actual binary search */
  for (limit = nmemb; limit != 0; limit >>= 1) {
    ptr = baseptr + (limit >> 1) * size;
    if ((rval = compar(key, ptr, data)) == 0) {
      /* element found */
      if (members)
        gt_array_add(members, ptr);
      else
        return ptr; /* because we got no <members> array */
      if (b) {
        /* mark found element */
        gt_bittab_set_bit(b, (ptr - (char *)base) / size);
      }
      /* looking left for equal elements */
      for (tmp_ptr = ptr - size;
           tmp_ptr >= (char *)base && !compar(key, tmp_ptr, data);
           tmp_ptr -= size) {
        gt_array_add(members, tmp_ptr);
        if (b)
          gt_bittab_set_bit(
            b, (tmp_ptr - (char *)base) / size); /* mark found element */
      }
      /* looking right for equal elements */
      for (tmp_ptr = ptr + size;
           tmp_ptr <(char *) base + nmemb * size && !compar(key, tmp_ptr, data);
           tmp_ptr += size) {
        gt_array_add(members, tmp_ptr);
        if (b)
          gt_bittab_set_bit(
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

void* gt_bsearch_data(const void *key, const void *base, size_t nmemb,
                      size_t size, GtCompareWithData compar, void *data)
{
  return bsearch_generic(NULL, key, base, nmemb, size, compar, data, NULL);
}

void gt_bsearch_all(GtArray *members, const void *key, const void *base,
                    size_t nmemb, size_t size, GtCompareWithData compar,
                    void *data)
{
  gt_assert(members);
  bsearch_generic(members, key, base, nmemb, size, compar, data, NULL);
}

void gt_bsearch_all_mark(GtArray *members, const void *key, const void *base,
                         size_t nmemb, size_t size, GtCompareWithData compar,
                         void *data, GtBittab *b)
{
  gt_assert(members);
  bsearch_generic(members, key, base, nmemb, size, compar, data, b);
}

static int cmp(const void *a_ptr, const void *b_ptr, GT_UNUSED void *unused)
{
  int a, b;
  gt_assert(a_ptr && b_ptr);
  a = *(int*) a_ptr;
  b = *(int*) b_ptr;
  if (a == b)
    return 0;
  if (a < b)
    return -1;
  return 1;
}

/* XXX: This unit test could be done much better by filling an array randomly,
   sorting it, and comparing gt_bsearch_all() against a brute force
   implementation.  */
int gt_bsearch_unit_test(GtError *err)
{
  GtArray *elements, *members;
  int key, element, *member_ptr;
  GtBittab *b;
  int had_err = 0;

  gt_error_check(err);

  elements = gt_array_new(sizeof (int));
  members = gt_array_new(sizeof (int*));

  /* the empty case */
  gt_bsearch_all(members, &key, gt_array_get_space(elements),
                 gt_array_size(elements), sizeof (int), cmp, NULL);
  gt_ensure(had_err, !gt_array_size(members)); /* no member found */
  gt_ensure(had_err, !gt_bsearch_data(&key, gt_array_get_space(elements),
                                   gt_array_size(elements), sizeof (int), cmp,
                                   NULL));

  /* 1 element */
  key = 7;
  element = 7;
  gt_array_add(elements, element);
  gt_bsearch_all(members, &key, gt_array_get_space(elements),
                 gt_array_size(elements), sizeof (int), cmp, NULL);
  gt_ensure(had_err, gt_array_size(members) == 1); /* one member found */
  member_ptr = *(int**) gt_array_get(members, 0);
  gt_ensure(had_err, *member_ptr == element);
  member_ptr = gt_bsearch_data(&key, gt_array_get_space(elements),
                               gt_array_size(elements), sizeof (int), cmp,
                               NULL);
  gt_ensure(had_err, *member_ptr == element);

  key = -7;
  gt_array_reset(members);
  gt_bsearch_all(members, &key, gt_array_get_space(elements),
                 gt_array_size(elements), sizeof (int), cmp, NULL);
  gt_ensure(had_err, !gt_array_size(members)); /* no member found */

  /* 2 elements */
  key = 7;
  gt_array_reset(members);
  gt_array_add(elements, element);
  gt_ensure(had_err, gt_array_size(elements) == 2);
  gt_bsearch_all(members, &key, gt_array_get_space(elements),
                 gt_array_size(elements), sizeof (int), cmp, NULL);
  gt_ensure(had_err, gt_array_size(members) == 2); /* two members found */
  member_ptr = *(int**) gt_array_get(members, 0);
  gt_ensure(had_err, *member_ptr == element);
  member_ptr = *(int**) gt_array_get(members, 1);
  gt_ensure(had_err, *member_ptr == element);
  gt_ensure(had_err, gt_bsearch_data(&key, gt_array_get_space(elements),
                                  gt_array_size(elements), sizeof (int), cmp,
                                  NULL));

  key = -7;
  gt_array_reset(members);
  gt_bsearch_all(members, &key, gt_array_get_space(elements),
                 gt_array_size(elements), sizeof (int), cmp, NULL);
  gt_ensure(had_err, !gt_array_size(members)); /* no member found */
  gt_ensure(had_err, !gt_bsearch_data(&key, gt_array_get_space(elements),
                                   gt_array_size(elements), sizeof (int), cmp,
                                   NULL));

  /* 3 elements */
  key = 7;
  gt_array_reset(members);
  gt_array_add(elements, element);
  gt_ensure(had_err, gt_array_size(elements) == 3);
  gt_bsearch_all(members, &key, gt_array_get_space(elements),
                 gt_array_size(elements), sizeof (int), cmp, NULL);
  gt_ensure(had_err, gt_array_size(members) == 3); /* three members found */
  member_ptr = *(int**) gt_array_get(members, 0);
  gt_ensure(had_err, *member_ptr == element);
  member_ptr = *(int**) gt_array_get(members, 1);
  gt_ensure(had_err, *member_ptr == element);
  member_ptr = *(int**) gt_array_get(members, 2);
  gt_ensure(had_err, *member_ptr == element);
  gt_ensure(had_err, gt_bsearch_data(&key, gt_array_get_space(elements),
                                  gt_array_size(elements), sizeof (int), cmp,
                                  NULL));

  key = -7;
  gt_array_reset(members);
  gt_bsearch_all(members, &key, gt_array_get_space(elements),
                 gt_array_size(elements), sizeof (int), cmp, NULL);
  gt_ensure(had_err, !gt_array_size(members)); /* no member found */
  gt_ensure(had_err, !gt_bsearch_data(&key, gt_array_get_space(elements),
                                   gt_array_size(elements), sizeof (int), cmp,
                                   NULL));

  /* large case: -10 -5 -3 -3 -3 0 1 2 3 */
  gt_array_reset(elements);
  element = -10;
  gt_array_add(elements, element);
  element = -5;
  gt_array_add(elements, element);
  element = -3;
  gt_array_add(elements, element);
  gt_array_add(elements, element);
  gt_array_add(elements, element);
  element = 0;
  gt_array_add(elements, element);
  element = 1;
  gt_array_add(elements, element);
  element = 2;
  gt_array_add(elements, element);
  element = 3;
  gt_array_add(elements, element);
  gt_ensure(had_err, gt_array_size(elements) == 9);
  key = -3;
  gt_array_reset(members);
  gt_bsearch_all(members, &key, gt_array_get_space(elements),
                 gt_array_size(elements), sizeof (int), cmp, NULL);
  gt_ensure(had_err, gt_array_size(members) == 3); /* three members found */
  member_ptr = *(int**) gt_array_get(members, 0);
  gt_ensure(had_err, *member_ptr == -3);
  member_ptr = *(int**) gt_array_get(members, 1);
  gt_ensure(had_err, *member_ptr == -3);
  member_ptr = *(int**) gt_array_get(members, 2);
  gt_ensure(had_err, *member_ptr == -3);
  gt_ensure(had_err,
         gt_bsearch_data(&key, gt_array_get_space(elements),
                         gt_array_size(elements), sizeof (int), cmp,NULL));

  /* test gt_bsearch_all_mark() with large case */
  gt_array_reset(members);
  b = gt_bittab_new(gt_array_size(elements));
  gt_bsearch_all_mark(members, &key, gt_array_get_space(elements),
                   gt_array_size(elements), sizeof (int), cmp, NULL, b);
  gt_ensure(had_err, gt_array_size(members) == 3); /* three members found */
  member_ptr = *(int**) gt_array_get(members, 0);
  gt_ensure(had_err, *member_ptr == -3);
  member_ptr = *(int**) gt_array_get(members, 1);
  gt_ensure(had_err, *member_ptr == -3);
  member_ptr = *(int**) gt_array_get(members, 2);
  gt_ensure(had_err, *member_ptr == -3);
  /* the correct elements are marked (and only these) */
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 2));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 3));
  gt_ensure(had_err, gt_bittab_bit_is_set(b, 4));
  gt_ensure(had_err, gt_bittab_count_set_bits(b) == 3);

  /* free */
  gt_array_delete(elements);
  gt_array_delete(members);
  gt_bittab_delete(b);

  return had_err;
}
