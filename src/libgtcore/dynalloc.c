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

#include <assert.h>
#include "libgtcore/dynalloc.h"
#include "libgtcore/ma.h"

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
  else {
    /* XXX: no overflow */
    assert(*allocated != SIZE_MAX);
    /* otherwise we double the allocated space, if possible */
    size_to_alloc = *allocated;
    while (size_to_alloc < size) {
      if (size_to_alloc > SIZE_MAX / 2) size_to_alloc = SIZE_MAX;
      else size_to_alloc *= 2;
    }
  }
  assert(size_to_alloc);
  rptr = ma_realloc(ptr, size_to_alloc);
  *allocated = size_to_alloc;
  return rptr;
}
