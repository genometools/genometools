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

#ifndef DLIST_H
#define DLIST_H

#include "libgtcore/error.h"
#include "libgtcore/fptr.h"

/* a double-linked list which is sorted according to a qsort(3)-like compare
   function (if one was supplied to the constructor) */
typedef struct Dlist Dlist;
typedef struct Dlistelem Dlistelem;

Dlist*        dlist_new(Compare);
Dlistelem*    dlist_first(const Dlist*);
Dlistelem*    dlist_last(const Dlist*);
Dlistelem*    dlist_find(const Dlist*, void*); /* O(n) */
unsigned long dlist_size(const Dlist*);
/* usually: O(n) (O(1) if data is added in sorted order) */
void          dlist_add(Dlist*, void *data);
/* frees the elem */
void          dlist_remove(Dlist*, Dlistelem*);
int           dlist_example(Error*);
int           dlist_unit_test(Error*);
void          dlist_delete(Dlist*);

Dlistelem*    dlistelem_next(const Dlistelem*);
Dlistelem*    dlistelem_previous(const Dlistelem*);
void*         dlistelem_get_data(const Dlistelem*);

#endif
