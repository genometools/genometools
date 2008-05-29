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

#ifndef COUNTINGSORT_H
#define COUNTINGSORT_H

#include <stdlib.h>
#include "libgtcore/error.h"

/*
  This module implements the counting sort algorithm. For a description see
  for example page 175 to page 177 of the book

  T.H. Cormen, C.E. Leiserson and R.L. Rivest. Introduction to Algorithms.
  MIT Press: Cambridge, MA, 1990.
*/

typedef unsigned long (*GetElemvalue)(const void *elem, void *data);

/* Sort the array of elements pointed to by <in> containing <size> many elements
   of size <elem_size> and store the result in the array <out> of the same size.
   <max_elemvalue> denotes the maximum value an element can have.
   <get_elemvalue> should return an integer value for the given element <elem>.
*/
void          countingsort(void *out, const void *in, size_t elem_size,
                           unsigned long size, unsigned long max_elemvalue,
                           void *data, GetElemvalue get_elemvalue);

/* If <max_elemvalue> is not known, it can be determined with this function. */
unsigned long countingsort_get_max(const void *in, size_t elem_size,
                                   unsigned long size, void *data,
                                   GetElemvalue get_elemvalue);

int           countingsort_unit_test(Error*);

#endif
