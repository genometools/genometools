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

#ifndef ARRAY_H
#define ARRAY_H

#include <stdlib.h>
#include "libgtcore/error.h"

typedef struct Array Array;

Array*        array_new(size_t);
Array*        array_clone(const Array*);
Array*        array_ref(Array*);
void*         array_get(const Array*, unsigned long);
void*         array_get_first(const Array*);
void*         array_get_last(const Array*);
void*         array_pop(Array*);
void*         array_get_space(const Array*);
#define       array_add(a, elem)\
              array_add_elem(a, &(elem), sizeof (elem))
void          array_add_elem(Array*, void*, size_t);
void          array_add_array(Array*, const Array*);
void          array_rem(Array*, unsigned long); /* O(n) */
void          array_reverse(Array*);
void          array_set_size(Array*, unsigned long);
void          array_reset(Array*);
size_t        array_elem_size(const Array*);
unsigned long array_size(const Array*);
void          array_sort(Array*, int(*compar)(const void*, const void*));
/* Compare the content of <array_a> with the content of <array_b>.
   <array_a> and <array_b> must have the same array_size() and
   array_elem_size(). */
int           array_cmp(const Array  *array_a, const Array *array_b);
int           array_iterate(const Array*,
                            int(*iterfunc)(void *info, const void *value,
                                           Error*),
                            void *info, Error*);
int           array_example(Error*);
int           array_unit_test(Error*);
void          array_delete(Array*);

#endif
