/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef COUNTINGSORT_H
#define COUNTINGSORT_H

#include <stdlib.h>
#include <libgt/env.h>

/*
  This module implements the counting sort algorithm. For a description see
  for example page 175 to page 177 of the book

  T.H. Cormen, C.E. Leiserson and R.L. Rivest. Introduction to Algorithms. MIT
  Press: Cambridge, MA, 1990.
*/

/* sort the array of elements pointed to by 'in' containing 'size' many elements
   of size 'elem_size' and store the result in the array 'out' of the same size.
   'max_elemvalue' denotes the maximum value an element can have.
   'get_elemvalue' should return an integer value for the given element 'elem'.
*/
void          countingsort(void *out, const void *in, size_t elem_size,
                           unsigned long size, unsigned long max_elemvalue,
                           void *data,
                           unsigned long (*get_elemvalue)(const void *elem,
                                                          void *data), Env*);

/* if 'max_elemvalue' is not known, it can be determined with this function */
unsigned long countingsort_get_max(const void *in, size_t elem_size,
                                   unsigned long size, void *data,
                                   unsigned long (*get_elemvalue)
                                                 (const void *elem,
                                                  void *data));

int           countingsort_unit_test(Env*);

#endif
