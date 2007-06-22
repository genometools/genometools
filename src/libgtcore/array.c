/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <limits.h>
#include <string.h>
#include <libgtcore/array.h>
#include <libgtcore/dynalloc.h>
#include <libgtcore/ensure.h>
#include <libgtcore/xansi.h>

#define NUM_OF_TESTS    100
#define MAX_SIZE        1024
#define UCCAST(PTR)             ((unsigned char *) (PTR))
#define MAKEPTR(PTR,ELEM,SIZE)  (UCCAST(PTR) + (ELEM) * (SIZE))

struct Array {
  void *space;
  unsigned long next_free;
  size_t allocated,
         size_of_elem;
};

Array* array_new(size_t size_of_elem, Env *env)
{
  Array *a = env_ma_calloc(env, 1, sizeof (Array));
  assert(size_of_elem);
  a->size_of_elem = size_of_elem;
  return a;
}

void* array_get(const Array *a, unsigned long idx)
{
  assert(a && idx < a->next_free);
  return MAKEPTR(a->space,idx,a->size_of_elem);
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
  return MAKEPTR(a->space,a->next_free,a->size_of_elem);
}

void array_rem(Array *a, unsigned long idx)
{
  unsigned long i;
  assert(a && idx < a->next_free);
  /* move elements */
  for (i = idx+1; i < a->next_free; i++) {
    memcpy(MAKEPTR(a->space,(i-1),a->size_of_elem),
           MAKEPTR(a->space,i,a->size_of_elem),
           a->size_of_elem);
  }
  /* remove last (now duplicated) element */
  a->next_free--;
}

void array_reverse(Array *a, Env *env)
{
  void *front, *back, *tmp;
  assert(a);
  tmp = env_ma_malloc(env, sizeof (a->size_of_elem));
  for (front = a->space, 
       back = MAKEPTR(a->space,a->next_free-1,a->size_of_elem);
       front < back; 
       UCCAST(front) += a->size_of_elem, UCCCAST(back) -= a->size_of_elem) {
    memcpy(tmp, front, a->size_of_elem);
    memcpy(front, back, a->size_of_elem);
    memcpy(back, tmp, a->size_of_elem);
  }
  env_ma_free(tmp, env);
}

void* array_get_space(const Array *a)
{
  assert(a);
  return a->space;
}

void array_add_elem(Array *a, void *elem, size_t size_of_elem, Env *env)
{
  assert(a && elem);
  assert(a->size_of_elem == size_of_elem);
  assert(a->next_free <= a->allocated);
  /* make sure we have enough space */
  if ((a->next_free + 1) * a->size_of_elem > a->allocated)
    a->space = dynalloc(a->space, &a->allocated,
                        (a->next_free + 1) * a->size_of_elem, env);
  /* add */
  memcpy(a->space + a->next_free * a->size_of_elem, elem, a->size_of_elem);
  a->next_free++;
}

void array_add_array(Array *dest, const Array *src, Env *env)
{
  unsigned long i;
  assert(dest && src && dest->size_of_elem == src->size_of_elem);
  for (i = 0; i < array_size(src); i++)
    array_add_elem(dest, array_get(src, i), src->size_of_elem, env);
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

void array_reset(Array *a)
{
  assert(a);
  a->next_free = 0;
}

Array* array_clone(const Array *a, Env *env)
{
  Array *a_copy;
  assert(a);
  a_copy = env_ma_malloc(env, sizeof (Array));
  /* XXX: overflow checks -> safemult(next_free, size_of_elem) */
  a_copy->space = env_ma_malloc(env, a->next_free * a->size_of_elem);
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
  char_array = array_new(sizeof (char), env);
  array_delete(char_array, env);
  int_array = array_new(sizeof (int), env);
  array_delete(int_array, env);

  char_array = array_new(sizeof (char), env);
  int_array = array_new(sizeof (int), env);
  char_array_test = env_ma_malloc(env, (MAX_SIZE + 1) * sizeof (char));
  int_array_test = env_ma_malloc(env, MAX_SIZE * sizeof (int));

  for (i = 0; i < NUM_OF_TESTS && !has_err; i++) {
    size = ((double) rand() / RAND_MAX) * MAX_SIZE;

    array_reset(char_array);
    array_reset(int_array);

    ensure(has_err, array_size(char_array) == 0);
    ensure(has_err, array_size(int_array) == 0);

    for (i = 0; i < size && !has_err; i++) {
      cc = ((double) rand() / RAND_MAX) * CHAR_MAX;
      ci = ((double) rand() / RAND_MAX) * INT_MAX;

      array_add(char_array, cc, env);
      array_add(int_array, ci, env);

      ensure(has_err, array_size(char_array) == i+1);
      ensure(has_err, array_size(int_array) == i+1);
      ensure(has_err, *((char*) array_get(char_array, i)) == cc);
      ensure(has_err, *((int*) array_get(int_array, i)) == ci);

      array_add_elem(char_array, &cc, sizeof (char), env);
      array_add_elem(int_array, &ci, sizeof (int), env);

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

  array_delete(char_array, env);
  array_delete(int_array, env);
  env_ma_free(char_array_test, env);
  env_ma_free(int_array_test, env);

  return has_err;
}

void array_delete(Array *a, Env *env)
{
  if (!a) return;
  env_ma_free(a->space, env);
  env_ma_free(a, env);
}
