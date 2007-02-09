/*
  Copyright (c) 2005-2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ARRAY_H
#define ARRAY_H

#include <stdlib.h>

typedef struct Array Array;

Array*        array_new(size_t);
Array*        array_clone(const Array*);
void*         array_get(const Array*, unsigned long);
void*         array_get_first(const Array*);
void*         array_get_last(const Array*);
void*         array_pop(Array*);
void*         array_get_space(const Array*);
#define       array_add(a, elem)\
              array_add_elem(a, &(elem), sizeof(elem))
void          array_add_elem(Array*, void*, size_t);
void          array_add_array(Array*, const Array*);
void          array_rem(Array*, unsigned long); /* O(n) */
void          array_reverse(Array*);
void          array_set_size(Array*, unsigned long);
size_t        array_elem_size(const Array*);
unsigned long array_size(const Array*);
int           array_unit_test(void);
void          array_free(Array*);

#if 0
  a typical array use case:

  int main(void)
  {
    unsigned long i;
    Array *a;
    a = array_new(sizeof(unsigned long));
    for (i = 0; i < 100; i++) {
      array_add(a, i);
      assert(i = *(unsigned long*) array_get(a, i));
    }
    assert(array_size(a) == 100);
    assert(*(unsigned long*) array_pop(a) == 99);
    assert(array_size(a) == 99);
    array_free(a);
    return EXIT_SUCCESS;
  }
#endif

#endif
