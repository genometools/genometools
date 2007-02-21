/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <limits.h>
#include "dlist.h"
#include "ensure.h"
#include "mathsupport.h"
#include "xansi.h"

#define NUM_OF_TESTS  100
#define MAX_SIZE      1024

struct Dlist {
  Compare cmp_func;
  Dlistelem *first,
            *last;
  unsigned long size;
};

struct Dlistelem {
  Dlistelem *previous,
            *next;
  void *data;
};

Dlist* dlist_new(Compare cmp_func)
{
  Dlist *dlist = xcalloc(1, sizeof (Dlist));
  dlist->cmp_func = cmp_func;
  return dlist;
}

Dlistelem* dlist_first(const Dlist *dlist)
{
  assert(dlist);
  return dlist->first;
}

Dlistelem* dlist_last(const Dlist *dlist)
{
  assert(dlist);
  return dlist->last;
}

unsigned long dlist_size(const Dlist *dlist)
{
  return dlist ? dlist->size : 0;
}

void dlist_add(Dlist *dlist, void *data)
{
  Dlistelem *oldelem, *newelem = xcalloc(1, sizeof (Dlistelem));
  assert(dlist); /* data can be null */
  newelem->data = data;

  if (!dlist->first) {
    assert(!dlist->last);
    dlist->first = newelem;
    dlist->last = newelem;
  }
  else {
    assert(dlist->first && dlist->last);
    if (dlist->cmp_func) {
      /* compare function defined -> find place to add new element */
      if (dlist->cmp_func(data, dlist->first->data) < 0) {
        /* the new element is smaller then the first element */
        assert(!dlist->first->previous);
        dlist->first->previous = newelem;
        newelem->next = dlist->first;
        dlist->first = newelem;
      }
      else if (dlist->cmp_func(dlist->last->data, data) <= 0) {
        /* the new element is larger or equal then the last element */
        assert(!dlist->last->next);
        dlist->last->next = newelem;
        newelem->previous = dlist->last;
        dlist->last = newelem;
      }
      else {
        /* traverse list backwards to find a place to insert the new element */
        oldelem = dlist->last->previous;
        assert(oldelem);
        while (oldelem) {
          if (dlist->cmp_func(oldelem->data, data) <= 0) {
            /* position found */
            assert(oldelem->next);
            newelem->next = oldelem->next;
            newelem->previous = oldelem;
            oldelem->next->previous = newelem;
            oldelem->next = newelem;
            break;
          }
          oldelem = oldelem->previous;
        }
        assert(oldelem); /* a position has been found */
      }
    }
    else {
      /* no compare function defined -> add new element to the end */
      assert(!dlist->last->next);
      dlist->last->next = newelem;
      newelem->previous = dlist->last;
      dlist->last = newelem;
    }
  }
  dlist->size++;
}

void dlist_remove(Dlist *dlist, Dlistelem *dlistelem)
{
  assert(dlist && dlistelem);
  assert(!dlistelem->previous || dlistelem->previous->next == dlistelem);
  assert(!dlistelem->next || dlistelem->next->previous == dlistelem);
  if (dlistelem->previous)
    dlistelem->previous->next = dlistelem->next;
  if (dlistelem->next)
    dlistelem->next->previous = dlistelem->previous;
  if (dlistelem == dlist->first)
    dlist->first = dlistelem->next;
  if (dlistelem == dlist->last)
    dlist->last = dlistelem->previous;
  dlist->size--;
  /* XXX: free(dlistelem); */
}

static int compare(const void *a, const void *b)
{
  return *(int*) a - *(int*) b;
}

int dlist_unit_test(Env *env)
{
  Dlist *dlist;
  Dlistelem *dlistelem;
  int i, j, size, *data,
      elem_a = 7,
      elem_b = 6,
      elems[MAX_SIZE],
      elems_backup[MAX_SIZE],
      has_err = 0;
  env_error_check(env);

  /* boundary case: empty dlist */
  dlist = dlist_new(compare);
  ensure(has_err, !dlist_size(dlist));
  dlist_delete(dlist);

  dlist = dlist_new(NULL);
  ensure(has_err, !dlist_size(dlist));
  dlist_delete(dlist);

  /* boundary case: dlist containing one element */
  dlist = dlist_new(compare);
  dlist_add(dlist, &elem_a);
  ensure(has_err, dlist_size(dlist) == 1);
  ensure(has_err, elem_a == *(int*) dlistelem_get_data(dlist_first(dlist)));
  dlist_delete(dlist);

  dlist = dlist_new(NULL);
  dlist_add(dlist, &elem_a);
  ensure(has_err, dlist_size(dlist) == 1);
  ensure(has_err, elem_a == *(int*) dlistelem_get_data(dlist_first(dlist)));
  dlist_delete(dlist);

  /* boundary case: dlist containing two elements */
  dlist = dlist_new(compare);
  dlist_add(dlist, &elem_a);
  dlist_add(dlist, &elem_b);
  ensure(has_err, dlist_size(dlist) == 2);
  ensure(has_err, elem_b == *(int*) dlistelem_get_data(dlist_first(dlist)));
  dlist_delete(dlist);

  dlist = dlist_new(NULL);
  dlist_add(dlist, &elem_a);
  dlist_add(dlist, &elem_b);
  ensure(has_err, dlist_size(dlist) == 2);
  ensure(has_err, elem_a == *(int*) dlistelem_get_data(dlist_first(dlist)));
  dlist_delete(dlist);

  for (i = 0; i < NUM_OF_TESTS && !has_err; i++) {
    /* construct the random elements for the list */
    size = rand_max(MAX_SIZE);
    for (j = 0; j < size; j++) {
      elems[j] = rand_max(INT_MAX);
      elems_backup[j] = elems[j];
    }

    /* sort the backup elements */
    qsort(elems_backup, size, sizeof (int), compare);

    /* test with compare function */
    dlist = dlist_new(compare);
    ensure(has_err, !dlist_size(dlist));
    for (j = 0; j < size && !has_err; j++) {
      dlist_add(dlist, elems + j);
      ensure(has_err, dlist_size(dlist) == j+1);

      for (dlistelem = dlist_first(dlist); dlistelem != NULL;
           dlistelem = dlistelem_next(dlistelem)) {
      }
    }
    j = 0;
    for (dlistelem = dlist_first(dlist); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      data = dlistelem_get_data(dlistelem);
      ensure(has_err, *data == elems_backup[j]);
      j++;
    }
    /* remove first element */
    if (dlist_size(dlist)) {
      dlist_remove(dlist, dlist_first(dlist));
      if (dlist_size(dlist)) {
        data = dlistelem_get_data(dlist_first(dlist));
        ensure(has_err, *data == elems_backup[1]);
      }
    }
    /* remove last element */
    if (dlist_size(dlist)) {
      dlist_remove(dlist, dlist_last(dlist));
      if (dlist_size(dlist)) {
        data = dlistelem_get_data(dlist_last(dlist));
        ensure(has_err, *data == elems_backup[size - 2]);
      }
    }
    /* XXX: fix this */
#if 0
    /* remove middle element */
    if (dlist_size(dlist) >= 2) {
      dlistelem = dlist_first(dlist);
      for (j = 1; j < dlist_size(dlist) / 2; j++)
        dlistelem = dlistelem_next(dlistelem);
      dlist_remove(dlist, dlistelem);
      dlistelem = dlist_first(dlist);
      for (j = 1; j < dlist_size(dlist) / 2 + 1; j++)
        dlistelem = dlistelem_next(dlistelem);
      data = dlistelem_get_data(dlist_last(dlist));
      ensure(has_err, *data == elems_backup[size / 2 + 1]);
    }
#endif
    dlist_delete(dlist);

    /* test without compare function */
    dlist = dlist_new(NULL);
    ensure(has_err, !dlist_size(dlist));
    for (j = 0; j < size && !has_err; j++) {
      dlist_add(dlist, elems + j);
      ensure(has_err, dlist_size(dlist) == j+1);
    }
    j = 0;
    for (dlistelem = dlist_first(dlist); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      data = dlistelem_get_data(dlistelem);
      ensure(has_err, *data == elems[j]);
      j++;
    }
    /* remove first element */
    if (dlist_size(dlist)) {
      dlist_remove(dlist, dlist_first(dlist));
      if (dlist_size(dlist)) {
        data = dlistelem_get_data(dlist_first(dlist));
        ensure(has_err, *data == elems[1]);
      }
    }
    /* remove last element */
    if (dlist_size(dlist)) {
      dlist_remove(dlist, dlist_last(dlist));
      if (dlist_size(dlist)) {
        data = dlistelem_get_data(dlist_last(dlist));
        ensure(has_err, *data == elems[size - 2]);
      }
    }
    dlist_delete(dlist);
  }

  return has_err;
}

void dlist_delete(Dlist *dlist)
{
  Dlistelem *elem;
  if (!dlist) return;
  elem = dlist->first;
  while (elem) {
    free(elem->previous);
    elem = elem->next;
  }
  free(dlist->last);
  free(dlist);
}

Dlistelem* dlistelem_next(const Dlistelem *dlistelem)
{
  assert(dlistelem);
  return dlistelem->next;
}

void* dlistelem_get_data(const Dlistelem *dlistelem)
{
  assert(dlistelem);
  return dlistelem->data;
}
