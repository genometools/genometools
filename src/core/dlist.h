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

#ifndef DLIST_H
#define DLIST_H

#include "core/error.h"
#include "core/fptr_api.h"

/* A double-linked list which is sorted according to a qsort(3)-like compare
   function (if one was supplied to the constructor). */
typedef struct GT_Dlist GT_Dlist;
typedef struct GT_Dlistelem GT_Dlistelem;

GT_Dlist*        gt_dlist_new(GT_Compare);
GT_Dlistelem*    gt_dlist_first(const GT_Dlist*);
GT_Dlistelem*    gt_dlist_last(const GT_Dlist*);
GT_Dlistelem*    gt_dlist_find(const GT_Dlist*, void*); /* O(n) */
unsigned long gt_dlist_size(const GT_Dlist*);
/* Usually O(n) (O(1) if data is added in sorted order). */
void          gt_dlist_add(GT_Dlist*, void *data);
/* Remove <dlistelem> from <dlist> and free it. */
void          gt_dlist_remove(GT_Dlist *dlist, GT_Dlistelem *dlistelem);
int           gt_dlist_example(GT_Error*);
int           gt_dlist_unit_test(GT_Error*);
void          gt_dlist_delete(GT_Dlist*);

GT_Dlistelem*    gt_dlistelem_next(const GT_Dlistelem*);
GT_Dlistelem*    gt_dlistelem_previous(const GT_Dlistelem*);
void*         gt_dlistelem_get_data(const GT_Dlistelem*);

#endif
