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

#include <assert.h>
#include <limits.h>
#include "libgtcore/dlist.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"

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
  Dlist *dlist = ma_calloc(1, sizeof (Dlist));
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

Dlistelem* dlist_find(const Dlist *dlist, void *new_data)
{
  Dlistelem *dlistelem;
  void *old_data;
  assert(dlist);
  for (dlistelem = dlist_first(dlist); dlistelem != NULL;
       dlistelem = dlistelem_next(dlistelem)) {
    old_data = dlistelem_get_data(dlistelem);
    if (dlist->cmp_func && !dlist->cmp_func(old_data, new_data))
      return dlistelem;
    else if (old_data == new_data)
      return dlistelem;
  }
  return NULL;
}

unsigned long dlist_size(const Dlist *dlist)
{
  return dlist ? dlist->size : 0;
}

void dlist_add(Dlist *dlist, void *data)
{
  Dlistelem *oldelem, *newelem;
  assert(dlist); /* data can be null */
  newelem = ma_calloc(1, sizeof (Dlistelem));
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
  ma_free(dlistelem);
}

static int intcompare(const void *a, const void *b)
{
  return *(int*) a - *(int*) b;
}

int dlist_example(UNUSED Error *err)
{
  Dlistelem *dlistelem;
  Dlist *dlist;
  void *data;
  int elem = 1984;
  error_check(err);

  dlist = dlist_new(NULL);
  dlist_add(dlist, &elem);
  dlist_add(dlist, &elem);
  dlist_add(dlist, &elem);

  /* a typical iterator loop */
  for (dlistelem = dlist_first(dlist); dlistelem != NULL;
       dlistelem = dlistelem_next(dlistelem)) {
    data = dlistelem_get_data(dlistelem);
    /* do something with data */
  }

  dlist_delete(dlist);

  return 0;
}

int dlist_unit_test(Error *err)
{
  Dlist *dlist;
  Dlistelem *dlistelem;
  int i, j, size, *data,
      elem_a = 7,
      elem_b = 6,
      elems[MAX_SIZE],
      elems_backup[MAX_SIZE],
      had_err = 0;
  error_check(err);

  /* boundary case: empty dlist */
  dlist = dlist_new(intcompare);
  ensure(had_err, !dlist_size(dlist));
  dlist_delete(dlist);

  dlist = dlist_new(NULL);
  ensure(had_err, !dlist_size(dlist));
  dlist_delete(dlist);

  /* boundary case: dlist containing one element */
  dlist = dlist_new(intcompare);
  dlist_add(dlist, &elem_a);
  ensure(had_err, dlist_size(dlist) == 1);
  ensure(had_err, elem_a == *(int*) dlistelem_get_data(dlist_first(dlist)));
  dlist_delete(dlist);

  dlist = dlist_new(NULL);
  dlist_add(dlist, &elem_a);
  ensure(had_err, dlist_size(dlist) == 1);
  ensure(had_err, elem_a == *(int*) dlistelem_get_data(dlist_first(dlist)));
  dlist_delete(dlist);

  /* boundary case: dlist containing two elements */
  dlist = dlist_new(intcompare);
  dlist_add(dlist, &elem_a);
  dlist_add(dlist, &elem_b);
  ensure(had_err, dlist_size(dlist) == 2);
  ensure(had_err, elem_b == *(int*) dlistelem_get_data(dlist_first(dlist)));
  dlist_delete(dlist);

  dlist = dlist_new(NULL);
  dlist_add(dlist, &elem_a);
  dlist_add(dlist, &elem_b);
  ensure(had_err, dlist_size(dlist) == 2);
  ensure(had_err, elem_a == *(int*) dlistelem_get_data(dlist_first(dlist)));
  dlist_delete(dlist);

  for (i = 0; i < NUM_OF_TESTS && !had_err; i++) {
    /* construct the random elements for the list */
    size = rand_max(MAX_SIZE);
    for (j = 0; j < size; j++) {
      elems[j] = rand_max(INT_MAX);
      elems_backup[j] = elems[j];
    }

    /* sort the backup elements */
    qsort(elems_backup, size, sizeof (int), intcompare);

    /* test with compare function */
    dlist = dlist_new(intcompare);
    ensure(had_err, !dlist_size(dlist));
    for (j = 0; j < size && !had_err; j++) {
      dlist_add(dlist, elems + j);
      ensure(had_err, dlist_size(dlist) == j+1);

      for (dlistelem = dlist_first(dlist); dlistelem != NULL;
           dlistelem = dlistelem_next(dlistelem)) {
      }
    }
    j = 0;
    for (dlistelem = dlist_first(dlist); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      data = dlistelem_get_data(dlistelem);
      ensure(had_err, *data == elems_backup[j]);
      j++;
    }
    /* test dlist_find() */
    for (j = 0; j < size; j++) {
      dlistelem = dlist_find(dlist, elems_backup + j);
      ensure(had_err, dlistelem);
      ensure(had_err, *(int*) dlistelem_get_data(dlistelem) == elems_backup[j]);
    }
    /* remove first element */
    if (dlist_size(dlist)) {
      dlist_remove(dlist, dlist_first(dlist));
      if (dlist_size(dlist)) {
        data = dlistelem_get_data(dlist_first(dlist));
        ensure(had_err, *data == elems_backup[1]);
      }
    }
    /* remove last element */
    if (dlist_size(dlist)) {
      dlist_remove(dlist, dlist_last(dlist));
      if (dlist_size(dlist)) {
        data = dlistelem_get_data(dlist_last(dlist));
        ensure(had_err, *data == elems_backup[size - 2]);
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
      ensure(had_err, *data == elems_backup[size / 2 + 1]);
    }
#endif
    dlist_delete(dlist);

    /* test without compare function */
    dlist = dlist_new(NULL);
    ensure(had_err, !dlist_size(dlist));
    for (j = 0; j < size && !had_err; j++) {
      dlist_add(dlist, elems + j);
      ensure(had_err, dlist_size(dlist) == j+1);
    }
    j = 0;
    for (dlistelem = dlist_first(dlist); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      data = dlistelem_get_data(dlistelem);
      ensure(had_err, *data == elems[j]);
      j++;
    }
    /* remove first element */
    if (dlist_size(dlist)) {
      dlist_remove(dlist, dlist_first(dlist));
      if (dlist_size(dlist)) {
        data = dlistelem_get_data(dlist_first(dlist));
        ensure(had_err, *data == elems[1]);
      }
    }
    /* remove last element */
    if (dlist_size(dlist)) {
      dlist_remove(dlist, dlist_last(dlist));
      if (dlist_size(dlist)) {
        data = dlistelem_get_data(dlist_last(dlist));
        ensure(had_err, *data == elems[size - 2]);
      }
    }
    dlist_delete(dlist);
  }

  return had_err;
}

void dlist_delete(Dlist *dlist)
{
  Dlistelem *elem;
  if (!dlist) return;
  elem = dlist->first;
  while (elem) {
    ma_free(elem->previous);
    elem = elem->next;
  }
  ma_free(dlist->last);
  ma_free(dlist);
}

Dlistelem* dlistelem_next(const Dlistelem *dlistelem)
{
  assert(dlistelem);
  return dlistelem->next;
}

Dlistelem* dlistelem_previous(const Dlistelem *dlistelem)
{
  assert(dlistelem);
  return dlistelem->previous;
}

void* dlistelem_get_data(const Dlistelem *dlistelem)
{
  assert(dlistelem);
  return dlistelem->data;
}
