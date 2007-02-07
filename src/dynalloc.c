/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "dynalloc.h"
#include "error.h"
#include "xansi.h"

#define SIZE_MAX ((size_t) ~0UL)

void* dynalloc(void *ptr, size_t *allocated, size_t size)
{
  size_t size_to_alloc = 0;
  void *rptr;
  assert(allocated && size);
  if (size <= *allocated)
    return ptr;
  if (*allocated == 0) {
    assert(ptr == NULL);
    /* if nothing has been allocated already, we allocate what was asked for */
    size_to_alloc = size;
  }
  else if (*allocated == SIZE_MAX) {
    /* we cannot allocate more space */
    error("overflow");
  }
  else {
    /* otherwise we double the allocated space, if possible */
    size_to_alloc = *allocated;
    while (size_to_alloc < size) {
      if (size_to_alloc > SIZE_MAX / 2) size_to_alloc = SIZE_MAX;
      else size_to_alloc *= 2;
    }
  }
  assert(size_to_alloc);
  rptr = xrealloc(ptr, size_to_alloc);
  *allocated = size_to_alloc;
  return rptr;
}
