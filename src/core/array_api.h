/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ARRAY_API_H
#define ARRAY_API_H

#include <stdlib.h>
#include "core/fptr_api.h"

typedef struct GT_Array GT_Array;

/* Return a new <GT_Array> whose elements have the size <size_of_elem>. */
GT_Array*     gt_array_new(size_t size_of_elem);
/* Return a clone of <array>. */
GT_Array*     gt_array_clone(const GT_Array *array);
/* Return pointer to element number <index> of <array>. <index> has to be
   smaller than <gt_array_size(array)>. */
void*         gt_array_get(const GT_Array *array, unsigned long index);
/* Return pointer to first element of <array>. */
void*         gt_array_get_first(const GT_Array *array);
/* Return pointer to last element of <array>. */
void*         gt_array_get_last(const GT_Array *array);
/* Return pointer to last element of <array> and remove it from <array>. */
void*         gt_array_pop(GT_Array *array);
/* Return pointer to the internal space of <array> where the elements are
   stored.  */
void*         gt_array_get_space(const GT_Array *array);
/* Add element <elem> to <array>. The size of <elem> must equal the given
   element size when the <array> was created and is determined automatically
   with the <sizeof> operator. */
#define       gt_array_add(array, elem) \
              gt_array_add_elem(array, &(elem), sizeof (elem))
/* Add element <elem> with size <size_of_elem> to <array>. <size_of_elem> must
   equal the given element size when the <array> was created. Usually, this
   method is not used directly and the macro <gt_array_add()> is used
   instead. */
void          gt_array_add_elem(GT_Array *array, void *elem,
                                size_t size_of_elem);
/* Add all elements of array <src> to the array <dest>. The element sizes of
   both arrays must be equal. */
void          gt_array_add_array(GT_Array *dest, const GT_Array *src);
/* Remove element with number <index> from <array> in O(<gt_array_size(array)>)
   time. <index> has to be smaller than <gt_array_size(array)>. */
void          gt_array_rem(GT_Array *array, unsigned long index);
/* Reverse the order of the elements in <array>. */
void          gt_array_reverse(GT_Array *array);
/* Set the size of <array> to <size>. <size> must be smaller than
   <gt_array_size(array)>. */
void          gt_array_set_size(GT_Array *array, unsigned long size);
/* Reset the <array>. That is, afterwards the array has size 0. */
void          gt_array_reset(GT_Array *array);
/* Return the size of the elements stored in <array>. */
size_t        gt_array_elem_size(const GT_Array *array);
/* Return the number of elements in <array>. If <array> equals <NULL>, 0 is
   returned. */
unsigned long gt_array_size(const GT_Array *array);
/* Sort <array> with the given compare function <compar>. */
void          gt_array_sort(GT_Array *array, GT_Compare compar);
/* Compare the content of <array_a> with the content of <array_b>.
   <array_a> and <array_b> must have the same gt_array_size() and
   gt_array_elem_size(). */
int           gt_array_cmp(const GT_Array *array_a, const GT_Array *array_b);
/* Delete <array>. */
void          gt_array_delete(GT_Array *array);

#endif
