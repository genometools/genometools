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

#ifndef DLIST_API_H
#define DLIST_API_H

#include "core/fptr_api.h"

/* A double-linked list which is sorted according to a <GtCompare> compare
   function (<qsort(3)>-like, only if one was supplied to the constructor). */
typedef struct GtDlist GtDlist;
typedef struct GtDlistelem GtDlistelem;

/* Creates a new <GtDlist> sorted according to the <GtCompare> function. If it
   is NULL, no sorting is enforced. */
GtDlist*      gt_dlist_new(GtCompare);

/* Returns the first <GtDlistelem> in a <GtDlist>. */
GtDlistelem*  gt_dlist_first(const GtDlist*);

/* Returns the last <GtDlistelem> in a <GtDlist>. */
GtDlistelem*  gt_dlist_last(const GtDlist*);

/* Returns the first <GtDlistelem> in a <GtDlist> which contains data identical
   to <data>. Takes O(n) time. */
GtDlistelem*  gt_dlist_find(const GtDlist*, void *data);

/* Returns the number of <GtDlistelem>s in a <GtDlist>. */
unsigned long gt_dlist_size(const GtDlist*);

/* Adds a new <GtDlistelem> containing <data> to a <GtDlist>. Usually O(n), but
   O(1) if data is added in sorted order. */
void          gt_dlist_add(GtDlist*, void *data);

/* Remove <dlistelem> from <dlist> and free it. */
void          gt_dlist_remove(GtDlist *dlist, GtDlistelem *dlistelem);

/* Example for usage of the <GtDlist> class. */
int           gt_dlist_example(GtError*);

/* Deletes a <GtDlist>. */
void          gt_dlist_delete(GtDlist*);

/* Returns the successor of a <GtDlistelem>, or NULL if the element is the last
   one in the <GtDlist>. */
GtDlistelem*  gt_dlistelem_next(const GtDlistelem*);

/* Returns the predecessor of a <GtDlistelem>, or NULL if the element is the
   first one in the <GtDlist>. */
GtDlistelem*  gt_dlistelem_previous(const GtDlistelem*);

/* Returns the data pointer attached to a <GtDlistelem>. */
void*         gt_dlistelem_get_data(const GtDlistelem*);

#endif
