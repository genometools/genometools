/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef BSEARCH_H
#define BSEARCH_H

#include "array.h"
#include "bittab.h"

/* similar interface to bsearch(3), except that all members which compare as
   equal are stored in the ``members'' array. The order in which the elements
   are added is undefined */
void bsearch_all(Array *members, const void *key, const void *base,
                 size_t nmemb, size_t size,
                 int (*compar) (const void *, const void *));

/* similar interface to bsearch_all(). Additionally, if a bittab is given (which
   must be of size ``nmemb''), the bits corresponding to the found elements are
   marked (i.e., set) */
void bsearch_all_mark(Array *members, const void *key, const void *base,
                      size_t nmemb, size_t size,
                      int (*compar) (const void *, const void *), Bittab*);

int bsearch_unit_test(Error*);

#endif
