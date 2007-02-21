/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <limits.h>
#include <string.h>
#include "array.h"
#include "dynalloc.h"
#include "ensure.h"
#include "xansi.h"

#define NUM_OF_TESTS	100
#define MAX_SIZE	1024

struct Array {
  void *space;
  unsigned long next_free;
  size_t allocated,
         size_of_elem;
};

Array* array_new(size_t size_of_elem)
{
  Array *a = xcalloc(1, sizeof (Array));
  assert(size_of_elem);
  a->size_of_elem = size_of_elem;
  return a;
}

void* array_get(const Array *a, unsigned long idx)
{
  assert(a && idx < a->next_free);
  return a->space + idx * a->size_of_elem;
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
  return a->space + a->next_free * a->size_of_elem;
}

void array_rem(Array *a, unsigned long idx)
{
  unsigned long i;
  assert(a && idx < a->next_free);
  /* move elements */
  for (i = idx+1; i < a->next_free; i++) {
    memcpy(a->space + (i-1) * a->size_of_elem, a->space + i * a->size_of_elem,
           a->size_of_elem);
  }
  /* remove last (now duplicated) element */
  a->next_free--;
}

void array_reverse(Array *a)
{
  void *front, *back, *tmp;
  assert(a);
  tmp = xmalloc(sizeof (a->size_of_elem));
  for (front = a->space, back = a->space + (a->next_free-1) * a->size_of_elem;
       front < back; front += a->size_of_elem, back -= a->size_of_elem) {
    memcpy(tmp, front, a->size_of_elem);
    memcpy(front, back, a->size_of_elem);
    memcpy(back, tmp, a->size_of_elem);
  }
  free(tmp);
}

void* array_get_space(const Array *a)
{
  assert(a);
  return a->space;
}

void array_add_elem(Array *a, void *elem, size_t size_of_elem)
{
  assert(a && elem);
  assert(a->size_of_elem == size_of_elem);
  assert(a->next_free <= a->allocated);
  /* make sure we have enough space */
  if ((a->next_free + 1) * a->size_of_elem > a->allocated)
    a->space = dynalloc(a->space, &a->allocated,
                        (a->next_free + 1) * a->size_of_elem);
  /* add */
  memcpy(a->space + a->next_free * a->size_of_elem, elem, a->size_of_elem);
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
  return a ? a->next_free :0;
}

void array_set_size(Array *a, unsigned long size)
{
  assert(a);
  assert(size <= a->next_free);
  a->next_free = size;
}

Array* array_clone(const Array *a)
{
  Array *a_copy;
  assert(a);
  a_copy = xmalloc(sizeof (Array));
  /* XXX: overflow checks -> safemult(next_free, size_of_elem) */
  a_copy->space = xmalloc(a->next_free * a->size_of_elem);
  memcpy(a_copy->space, a->space, a->next_free * a->size_of_elem);
  a_copy->next_free = a_copy->allocated = a->next_free;
  a_copy->size_of_elem = a->size_of_elem;
  return a_copy;
}

int array_unit_test(Env *env)
{
  Array *char_array, *int_array;
  char cc, *char_array_test;
  int ci, *int_array_test;
  unsigned long i, j, size;
  int has_err = 0;
  env_error_check(env);

  /* testing an empty array */
  char_array = array_new(sizeof (char));
  array_delete(char_array);
  int_array = array_new(sizeof (int));
  array_delete(int_array);

  char_array = array_new(sizeof (char));
  int_array = array_new(sizeof (int));
  char_array_test = xmalloc((MAX_SIZE + 1) * sizeof (char));
  int_array_test = xmalloc(MAX_SIZE * sizeof (int));

  for (i = 0; i < NUM_OF_TESTS && !has_err; i++) {
    size = ((double) rand() / RAND_MAX) * MAX_SIZE;

    array_set_size(char_array, 0);
    array_set_size(int_array, 0);

    ensure(has_err, array_size(char_array) == 0);
    ensure(has_err, array_size(int_array) == 0);

    for (i = 0; i < size && !has_err; i++) {
      cc = ((double) rand() / RAND_MAX) * CHAR_MAX;
      ci = ((double) rand() / RAND_MAX) * INT_MAX;

      array_add(char_array, cc);
      array_add(int_array, ci);

      ensure(has_err, array_size(char_array) == i+1);
      ensure(has_err, array_size(int_array) == i+1);
      ensure(has_err, *((char*) array_get(char_array, i)) == cc);
      ensure(has_err, *((int*) array_get(int_array, i)) == ci);

      array_add_elem(char_array, &cc, sizeof (char));
      array_add_elem(int_array, &ci, sizeof (int));

      ensure(has_err, array_size(char_array) == i+2);
      ensure(has_err, array_size(int_array) == i+2);
      ensure(has_err, *((char*) array_get(char_array, i+1)) == cc);
      ensure(has_err, *((int*) array_get(int_array, i+1)) == ci);
      ensure(has_err, *((char*) array_pop(char_array)) == cc);
      ensure(has_err, *((int*) array_pop(int_array)) == ci);
      ensure(has_err, array_size(char_array) == i+1);
      ensure(has_err, array_size(int_array) == i+1);
      ensure(has_err, *((char*) array_get(char_array, i)) == cc);
      ensure(has_err, *((int*) array_get(int_array, i)) == ci);

      char_array_test[i] = cc;
      char_array_test[i+1]= '\0';
      int_array_test[i] = ci;

      ensure(has_err, strncmp(array_get_space(char_array), char_array_test,
                              strlen(char_array_test)) == 0);

      for (j = 0; j <= i && !has_err; j++)
        ensure(has_err, *(int*) array_get(int_array, j) == int_array_test[j]);
    }
  }

  array_delete(char_array);
  array_delete(int_array);
  free(char_array_test);
  free(int_array_test);

  return has_err;
}

void array_delete(Array *a)
{
  if (!a) return;
  free(a->space);
  free(a);
}
