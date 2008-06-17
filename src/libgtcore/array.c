/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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
#include <limits.h>
#include <string.h>
#include "libgtcore/array.h"
#include "libgtcore/dynalloc.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtcore/range.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"

#define NUM_OF_TESTS    100
#define MAX_SIZE        1024

struct Array {
  void *space;
  unsigned long next_free;
  size_t allocated,
         size_of_elem;
  unsigned int reference_count;
};

Array* array_new(size_t size_of_elem)
{
  Array *a = ma_calloc(1, sizeof (Array));
  assert(size_of_elem);
  a->size_of_elem = size_of_elem;
  return a;
}

Array* array_clone(const Array *a)
{
  Array *a_copy;
  assert(a);
  a_copy = ma_malloc(sizeof (Array));
  a_copy->space = ma_malloc(a->next_free * a->size_of_elem);
  memcpy(a_copy->space, a->space, a->next_free * a->size_of_elem);
  a_copy->next_free = a_copy->allocated = a->next_free;
  a_copy->size_of_elem = a->size_of_elem;
  a_copy->reference_count = 0;
  return a_copy;
}

Array* array_ref(Array *a)
{
  if (!a) return NULL;
  a->reference_count++;
  return a;
}

void* array_get(const Array *a, unsigned long idx)
{
  assert(a && idx < a->next_free);
  return (char*) a->space + idx * a->size_of_elem;
}

void* array_get_first(const Array *a)
{
  return array_get(a, 0);
}

void* array_get_last(const Array *a)
{
  assert(a->next_free);
  return array_get(a, a->next_free-1);
}

void* array_pop(Array *a)
{
  assert(a && a->next_free);
  a->next_free--;
  return (char*) a->space + a->next_free * a->size_of_elem;
}

void array_rem(Array *a, unsigned long idx)
{
  unsigned long i;
  assert(a && idx < a->next_free);
  /* move elements */
  for (i = idx+1; i < a->next_free; i++) {
    memcpy((char*) a->space + (i-1) * a->size_of_elem,
           (char*) a->space + i * a->size_of_elem,
           a->size_of_elem);
  }
  /* remove last (now duplicated) element */
  a->next_free--;
}

void array_reverse(Array *a)
{
  char *front, *back, *tmp;
  assert(a);
  tmp = ma_malloc(a->size_of_elem);
  for (front = a->space,
       back = (char*) a->space + (a->next_free-1) * a->size_of_elem;
       front < back;
       front += a->size_of_elem, back -= a->size_of_elem) {
    memcpy(tmp, front, a->size_of_elem);
    memcpy(front, back, a->size_of_elem);
    memcpy(back, tmp, a->size_of_elem);
  }
  ma_free(tmp);
}

void* array_get_space(const Array *a)
{
  assert(a);
  return a->space;
}

void array_add_elem(Array *a, void *elem, UNUSED size_t size_of_elem)
{
  assert(a && elem);
  assert(a->size_of_elem == size_of_elem);
  assert(a->next_free <= a->allocated);
  /* make sure we have enough space */
  if ((a->next_free + 1) * a->size_of_elem > a->allocated) {
    a->space = dynalloc(a->space, &a->allocated,
                        (a->next_free + 1) * a->size_of_elem);
  }
  /* add */
  memcpy((char*) a->space + a->next_free * a->size_of_elem, elem,
         a->size_of_elem);
  a->next_free++;
}

void array_add_array(Array *dest, const Array *src)
{
  unsigned long i;
  assert(dest && src && dest->size_of_elem == src->size_of_elem);
  for (i = 0; i < array_size(src); i++)
    array_add_elem(dest, array_get(src, i), src->size_of_elem);
}

size_t array_elem_size(const Array *a)
{
  assert(a);
  return a->size_of_elem;
}

unsigned long array_size(const Array *a)
{
  return a ? a->next_free : 0;
}

void array_set_size(Array *a, unsigned long size)
{
  assert(a);
  assert(size <= a->next_free);
  a->next_free = size;
}

void array_reset(Array *a)
{
  assert(a);
  a->next_free = 0;
}

void array_sort(Array *a, int(*compar)(const void*, const void*))
{
  assert(a && compar);
  qsort(a->space, a->next_free, a->size_of_elem, compar);
}

int array_cmp(const Array *array_a, const Array *array_b)
{
  assert(array_size(array_a) == array_size(array_b));
  assert(array_elem_size(array_a) == array_elem_size(array_b));
  return memcmp(array_a->space, array_b->space,
                array_a->size_of_elem * array_a->next_free);
}

int array_iterate(const Array *a,
                  int(*iterfunc)(void *info, const void *value, Error *err),
                  void *info, Error *err)
{
  unsigned long idx;
  error_check(err);
  assert(a && iterfunc);
  for (idx = 0; idx < array_size(a); idx++) {
    if (iterfunc(info, array_get(a, idx), err))
      return -1;
  }
  return 0;
}

static int iterate_test_func(void *info, const void *value, UNUSED Error *err)
{
  unsigned long *i;
  Range range;
  int had_err = 0;
  error_check(err);
  i = (unsigned long*) info;
  range = *(Range*) value;
  /* XXX: uncomment */
#if 0
  ensure(had_err, range.start == *i + 1);
  ensure(had_err, range.end == *i + 101);
#endif
  (*i)++;
  return had_err;
}

static int iterate_fail_func(UNUSED void *info, UNUSED const void *value,
                             UNUSED Error *err)
{
  return -1;
}

int array_example(UNUSED Error *err)
{
  unsigned long i;
  Array *a;

  error_check(err);

  /* an example array use case */

  a = array_new(sizeof (unsigned long));
  for (i = 0; i < 100; i++) {
    array_add(a, i);
    assert(i == *(unsigned long*) array_get(a, i));
  }
  assert(array_size(a) == 100);
  assert(*(unsigned long*) array_pop(a) == 99);
  assert(array_size(a) == 99);
  array_delete(a);

  return 0;
}

int array_unit_test(Error *err)
{
  Array *char_array, *int_array, *a = NULL, *aref, *aclone;
  char cc, *char_array_test;
  int ci, *int_array_test;
  unsigned long i, j, size;
  Range range;
  int had_err = 0;
  error_check(err);

  /* testing an empty array */
  char_array = array_new(sizeof (char));
  array_delete(char_array);
  int_array = array_new(sizeof (int));
  array_delete(int_array);

  char_array = array_new(sizeof (char));
  int_array = array_new(sizeof (int));
  char_array_test = ma_malloc((MAX_SIZE + 1) * sizeof (char));
  int_array_test = ma_malloc(MAX_SIZE * sizeof (int));

  for (i = 0; !had_err && i < NUM_OF_TESTS; i++) {
    size = rand_max(MAX_SIZE);

    array_reset(char_array);
    array_set_size(int_array, 0);

    ensure(had_err, array_size(char_array) == 0);
    ensure(had_err, array_size(int_array) == 0);

    for (i = 0; !had_err && i < size; i++) {
      cc = rand_max(CHAR_MAX);
      ci = rand_max(INT_MAX);

      array_add(char_array, cc);
      array_add(int_array, ci);

      ensure(had_err, array_size(char_array) == i+1);
      ensure(had_err, array_size(int_array) == i+1);
      ensure(had_err, *((char*) array_get(char_array, i)) == cc);
      ensure(had_err, *((int*) array_get(int_array, i)) == ci);

      array_add_elem(char_array, &cc, sizeof (char));
      array_add_elem(int_array, &ci, sizeof (int));

      ensure(had_err, array_size(char_array) == i+2);
      ensure(had_err, array_size(int_array) == i+2);
      ensure(had_err, *((char*) array_get(char_array, i+1)) == cc);
      ensure(had_err, *((int*) array_get(int_array, i+1)) == ci);
      ensure(had_err, *((char*) array_pop(char_array)) == cc);
      ensure(had_err, *((int*) array_pop(int_array)) == ci);
      ensure(had_err, array_size(char_array) == i+1);
      ensure(had_err, array_size(int_array) == i+1);
      ensure(had_err, *((char*) array_get(char_array, i)) == cc);
      ensure(had_err, *((int*) array_get(int_array, i)) == ci);

      char_array_test[i] = cc;
      char_array_test[i+1]= '\0';
      int_array_test[i] = ci;

      ensure(had_err, strncmp(array_get_space(char_array), char_array_test,
                              strlen(char_array_test)) == 0);

      for (j = 0; j <= i && !had_err; j++)
        ensure(had_err, *(int*) array_get(int_array, j) == int_array_test[j]);
    }
  }

  /* test array_reverse(), array_iterate(), array_rem(), and array_ref() */
  if (!had_err) {
    a = array_new(sizeof (Range));
    for (i = 0; i < 24; i++) {
      range.start = i + 1;
      range.end   = i + 101;
      array_add(a, range);
    }
    array_reverse(a);
    for (i = 0; !had_err && i < 24; i++) {
      range.start = i + 1;
      range.end   = i + 101;
      ensure(had_err, !range_compare(range, *(Range*) array_get(a, 23 - i)));
    }
    aref = array_ref(a);
    array_delete(aref);
  }
  if (!had_err) {
    array_reverse(a);
    i = 0;
    ensure(had_err, !array_iterate(a, iterate_test_func, &i, err));
    ensure(had_err, array_iterate(a, iterate_fail_func, NULL, err));
  }
  if (!had_err) {
    aclone = array_clone(a);
    for (i = 0;!had_err && i < array_size(a); i++) {
      ensure(had_err, !range_compare_ptr(array_get(a, i),
                                         array_get(aclone, i)));
    }
    array_delete(aclone);
  }
  if (!had_err) {
    array_rem(a, 13);
    array_rem(a, 12);
    array_rem(a, 11);
    array_rem(a, 10);
    ensure(had_err, array_size(a) == 20);
    for (i = 0; !had_err && i < 10; i++) {
      range.start = i + 1;
      range.end   = i + 101;
      ensure(had_err, !range_compare(range, *(Range*) array_get(a, i)));
    }
    for (i = 10; !had_err && i < 20; i++) {
      range.start = 4 + i + 1;
      range.end   = 4 + i + 101;
      ensure(had_err, !range_compare(range, *(Range*) array_get(a, i)));
    }
  }
  array_delete(a);

  array_delete(char_array);
  array_delete(int_array);
  ma_free(char_array_test);
  ma_free(int_array_test);

  return had_err;
}

void array_delete(Array *a)
{
  if (!a) return;
  if (a->reference_count) {
    a->reference_count--;
    return;
  }
  ma_free(a->space);
  ma_free(a);
}
