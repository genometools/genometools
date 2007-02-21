/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef DLIST_H
#define DLIST_H

#include "env.h"
#include "fptr.h"

/* a double-linked list which is sorted according to a qsort(3)-like compare
   function (if one was supplied to the constructor) */
typedef struct Dlist Dlist;
typedef struct Dlistelem Dlistelem;

Dlist*        dlist_new(Compare);
Dlistelem*    dlist_first(const Dlist*);
Dlistelem*    dlist_last(const Dlist*);
unsigned long dlist_size(const Dlist*);
void          dlist_add(Dlist*, void *data); /* usually: O(n) (O(1) if data is
                                                added in sorted order) */
void          dlist_remove(Dlist*, Dlistelem*); /* XXX: frees the elem */
int           dlist_unit_test(Env*);
void          dlist_delete(Dlist*);

Dlistelem*    dlistelem_next(const Dlistelem*);
void*         dlistelem_get_data(const Dlistelem*);

#if 0
  a typical iterator loop:

  for (dlistelem = dlist_first(dlist); dlistelem != NULL;
       dlistelem = dlistelem_next(dlistelem)) {
    data = dlistelem_get_data(dlistelem);
    /* do something with data */
  }
#endif

#endif
