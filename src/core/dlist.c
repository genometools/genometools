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

#include <limits.h>
#include "core/dlist.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"

#define NUM_OF_TESTS  100
#define MAX_SIZE      1024

struct GtDlist {
  GtCompareWithData cmp_func;
  GtDlistelem *first,
              *last;
  void *data;
  unsigned long size;
};

struct GtDlistelem {
  GtDlistelem *previous,
              *next;
  void *data;
};

static int gt_dlist_cmp_wrapper(const void *a, const void *b, void *data)
{
  return ((GtCompare) data)(a, b);
}

GtDlist* gt_dlist_new(GtCompare cmp_func)
{
  GtDlist *dlist = gt_calloc(1, sizeof (GtDlist));
  if (cmp_func == NULL)
    dlist->cmp_func = NULL;
  else
    dlist->cmp_func = gt_dlist_cmp_wrapper;
  dlist->data = cmp_func;
  return dlist;
}

GtDlist* gt_dlist_new_with_data(GtCompareWithData cmp_func, void *data)
{
  GtDlist *dlist = gt_calloc(1, sizeof (GtDlist));
  dlist->cmp_func = cmp_func;
  dlist->data = data;
  return dlist;
}

GtDlistelem* gt_dlist_first(const GtDlist *dlist)
{
  gt_assert(dlist);
  return dlist->first;
}

GtDlistelem* gt_dlist_last(const GtDlist *dlist)
{
  gt_assert(dlist);
  return dlist->last;
}

GtDlistelem* gt_dlist_find(const GtDlist *dlist, void *new_data)
{
  GtDlistelem *dlistelem;
  void *old_data;
  gt_assert(dlist);
  for (dlistelem = gt_dlist_first(dlist); dlistelem != NULL;
       dlistelem = gt_dlistelem_next(dlistelem)) {
    old_data = gt_dlistelem_get_data(dlistelem);
    if (dlist->cmp_func && !dlist->cmp_func(old_data, new_data, dlist->data))
      return dlistelem;
    else if (old_data == new_data)
      return dlistelem;
  }
  return NULL;
}

unsigned long gt_dlist_size(const GtDlist *dlist)
{
  return dlist ? dlist->size : 0;
}

void gt_dlist_add(GtDlist *dlist, void *data)
{
  GtDlistelem *oldelem, *newelem;
  gt_assert(dlist); /* data can be null */
  newelem = gt_calloc(1, sizeof (GtDlistelem));
  newelem->data = data;

  if (!dlist->first) {
    gt_assert(!dlist->last);
    dlist->first = newelem;
    dlist->last = newelem;
  }
  else {
    gt_assert(dlist->first && dlist->last);
    if (dlist->cmp_func) {
      /* compare function defined -> find place to add new element */
      if (dlist->cmp_func(data, dlist->first->data, dlist->data) < 0) {
        /* the new element is smaller then the first element */
        gt_assert(!dlist->first->previous);
        dlist->first->previous = newelem;
        newelem->next = dlist->first;
        dlist->first = newelem;
      }
      else if (dlist->cmp_func(dlist->last->data, data, dlist->data) <= 0) {
        /* the new element is larger or equal then the last element */
        gt_assert(!dlist->last->next);
        dlist->last->next = newelem;
        newelem->previous = dlist->last;
        dlist->last = newelem;
      }
      else {
        /* traverse list backwards to find a place to insert the new element */
        oldelem = dlist->last->previous;
        gt_assert(oldelem);
        while (oldelem) {
          if (dlist->cmp_func(oldelem->data, data, dlist->data) <= 0) {
            /* position found */
            gt_assert(oldelem->next);
            newelem->next = oldelem->next;
            newelem->previous = oldelem;
            oldelem->next->previous = newelem;
            oldelem->next = newelem;
            break;
          }
          oldelem = oldelem->previous;
        }
        gt_assert(oldelem); /* a position has been found */
      }
    }
    else {
      /* no compare function defined -> add new element to the end */
      gt_assert(!dlist->last->next);
      dlist->last->next = newelem;
      newelem->previous = dlist->last;
      dlist->last = newelem;
    }
  }
  dlist->size++;
}

void gt_dlist_remove(GtDlist *dlist, GtDlistelem *dlistelem)
{
  gt_assert(dlist && dlistelem);
  gt_assert(!dlistelem->previous || dlistelem->previous->next == dlistelem);
  gt_assert(!dlistelem->next || dlistelem->next->previous == dlistelem);
  if (dlistelem->previous)
    dlistelem->previous->next = dlistelem->next;
  if (dlistelem->next)
    dlistelem->next->previous = dlistelem->previous;
  if (dlistelem == dlist->first)
    dlist->first = dlistelem->next;
  if (dlistelem == dlist->last)
    dlist->last = dlistelem->previous;
  dlist->size--;
  gt_free(dlistelem);
}

static int intcompare(const void *a, const void *b)
{
  return *(int*) a - *(int*) b;
}

int gt_dlist_example(GT_UNUSED GtError *err)
{
  GtDlistelem *dlistelem;
  GtDlist *dlist;
  GT_UNUSED void *data;
  int elem = 1984;
  gt_error_check(err);

  dlist = gt_dlist_new(NULL);
  gt_dlist_add(dlist, &elem);
  gt_dlist_add(dlist, &elem);
  gt_dlist_add(dlist, &elem);

  /* a typical iterator loop */
  for (dlistelem = gt_dlist_first(dlist); dlistelem != NULL;
       dlistelem = gt_dlistelem_next(dlistelem)) {
    data = gt_dlistelem_get_data(dlistelem);
    /* do something with data */
  }

  gt_dlist_delete(dlist);

  return 0;
}

int gt_dlist_unit_test(GtError *err)
{
  GtDlist *dlist;
  GtDlistelem *dlistelem;
  int i, j, size, *data,
      elem_a = 7,
      elem_b = 6,
      elems[MAX_SIZE],
      elems_backup[MAX_SIZE],
      had_err = 0;
  gt_error_check(err);

  /* boundary case: empty dlist */
  dlist = gt_dlist_new(intcompare);
  gt_ensure(had_err, !gt_dlist_size(dlist));
  gt_dlist_delete(dlist);

  dlist = gt_dlist_new(NULL);
  gt_ensure(had_err, !gt_dlist_size(dlist));
  gt_dlist_delete(dlist);

  /* boundary case: dlist containing one element */
  dlist = gt_dlist_new(intcompare);
  gt_dlist_add(dlist, &elem_a);
  gt_ensure(had_err, gt_dlist_size(dlist) == 1);
  gt_ensure(had_err,
         elem_a == *(int*) gt_dlistelem_get_data(gt_dlist_first(dlist)));
  gt_dlist_delete(dlist);

  dlist = gt_dlist_new(NULL);
  gt_dlist_add(dlist, &elem_a);
  gt_ensure(had_err, gt_dlist_size(dlist) == 1);
  gt_ensure(had_err,
         elem_a == *(int*) gt_dlistelem_get_data(gt_dlist_first(dlist)));
  gt_dlist_delete(dlist);

  /* boundary case: dlist containing two elements */
  dlist = gt_dlist_new(intcompare);
  gt_dlist_add(dlist, &elem_a);
  gt_dlist_add(dlist, &elem_b);
  gt_ensure(had_err, gt_dlist_size(dlist) == 2);
  gt_ensure(had_err,
         elem_b == *(int*) gt_dlistelem_get_data(gt_dlist_first(dlist)));
  gt_dlist_delete(dlist);

  dlist = gt_dlist_new(NULL);
  gt_dlist_add(dlist, &elem_a);
  gt_dlist_add(dlist, &elem_b);
  gt_ensure(had_err, gt_dlist_size(dlist) == 2);
  gt_ensure(had_err,
         elem_a == *(int*) gt_dlistelem_get_data(gt_dlist_first(dlist)));
  gt_dlist_delete(dlist);

  for (i = 0; i < NUM_OF_TESTS && !had_err; i++) {
    /* construct the random elements for the list */
    size = gt_rand_max(MAX_SIZE);
    for (j = 0; j < size; j++) {
      elems[j] = gt_rand_max(INT_MAX);
      elems_backup[j] = elems[j];
    }

    /* sort the backup elements */
    qsort(elems_backup, size, sizeof (int), intcompare);

    /* test with compare function */
    dlist = gt_dlist_new(intcompare);
    gt_ensure(had_err, !gt_dlist_size(dlist));
    for (j = 0; j < size && !had_err; j++) {
      gt_dlist_add(dlist, elems + j);
      gt_ensure(had_err, gt_dlist_size(dlist) == j+1);

      for (dlistelem = gt_dlist_first(dlist); dlistelem != NULL;
           dlistelem = gt_dlistelem_next(dlistelem)) {
      }
    }
    j = 0;
    for (dlistelem = gt_dlist_first(dlist); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      data = gt_dlistelem_get_data(dlistelem);
      gt_ensure(had_err, *data == elems_backup[j]);
      j++;
    }
    /* test gt_dlist_find() */
    for (j = 0; j < size; j++) {
      dlistelem = gt_dlist_find(dlist, elems_backup + j);
      gt_ensure(had_err, dlistelem);
      gt_ensure(had_err,
             *(int*) gt_dlistelem_get_data(dlistelem) == elems_backup[j]);
    }
    /* remove first element */
    if (gt_dlist_size(dlist)) {
      gt_dlist_remove(dlist, gt_dlist_first(dlist));
      if (gt_dlist_size(dlist)) {
        data = gt_dlistelem_get_data(gt_dlist_first(dlist));
        gt_ensure(had_err, *data == elems_backup[1]);
      }
    }
    /* remove last element */
    if (gt_dlist_size(dlist)) {
      gt_dlist_remove(dlist, gt_dlist_last(dlist));
      if (gt_dlist_size(dlist)) {
        data = gt_dlistelem_get_data(gt_dlist_last(dlist));
        gt_ensure(had_err, *data == elems_backup[size - 2]);
      }
    }
    /* XXX: fix this */
#if 0
    /* remove middle element */
    if (gt_dlist_size(dlist) >= 2) {
      dlistelem = gt_dlist_first(dlist);
      for (j = 1; j < gt_dlist_size(dlist) / 2; j++)
        dlistelem = gt_dlistelem_next(dlistelem);
      gt_dlist_remove(dlist, dlistelem);
      dlistelem = gt_dlist_first(dlist);
      for (j = 1; j < gt_dlist_size(dlist) / 2 + 1; j++)
        dlistelem = gt_dlistelem_next(dlistelem);
      data = gt_dlistelem_get_data(gt_dlist_last(dlist));
      gt_ensure(had_err, *data == elems_backup[size / 2 + 1]);
    }
#endif
    gt_dlist_delete(dlist);

    /* test without compare function */
    dlist = gt_dlist_new(NULL);
    gt_ensure(had_err, !gt_dlist_size(dlist));
    for (j = 0; j < size && !had_err; j++) {
      gt_dlist_add(dlist, elems + j);
      gt_ensure(had_err, gt_dlist_size(dlist) == j+1);
    }
    j = 0;
    for (dlistelem = gt_dlist_first(dlist); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      data = gt_dlistelem_get_data(dlistelem);
      gt_ensure(had_err, *data == elems[j]);
      j++;
    }
    /* remove first element */
    if (gt_dlist_size(dlist)) {
      gt_dlist_remove(dlist, gt_dlist_first(dlist));
      if (gt_dlist_size(dlist)) {
        data = gt_dlistelem_get_data(gt_dlist_first(dlist));
        gt_ensure(had_err, *data == elems[1]);
      }
    }
    /* remove last element */
    if (gt_dlist_size(dlist)) {
      gt_dlist_remove(dlist, gt_dlist_last(dlist));
      if (gt_dlist_size(dlist)) {
        data = gt_dlistelem_get_data(gt_dlist_last(dlist));
        gt_ensure(had_err, *data == elems[size - 2]);
      }
    }
    gt_dlist_delete(dlist);
  }

  return had_err;
}

void gt_dlist_delete(GtDlist *dlist)
{
  GtDlistelem *elem;
  if (!dlist) return;
  elem = dlist->first;
  while (elem) {
    gt_free(elem->previous);
    elem = elem->next;
  }
  gt_free(dlist->last);
  gt_free(dlist);
}

GtDlistelem* gt_dlistelem_next(const GtDlistelem *dlistelem)
{
  gt_assert(dlistelem);
  return dlistelem->next;
}

GtDlistelem* gt_dlistelem_previous(const GtDlistelem *dlistelem)
{
  gt_assert(dlistelem);
  return dlistelem->previous;
}

void* gt_dlistelem_get_data(const GtDlistelem *dlistelem)
{
  gt_assert(dlistelem);
  return dlistelem->data;
}
