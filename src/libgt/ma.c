/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdbool.h>
#include "ma.h"
#include "xansi.h"

/* the memory allocator class */
struct MA {
  bool bookkeeping;
};

MA* ma_new(void)
{
  return xcalloc(1, sizeof (MA));
}

void* ma_malloc(MA *ma, size_t size)
{
  assert(ma);
  return xmalloc(size);
}

void* ma_calloc(MA *ma, size_t nmemb, size_t size)
{
  assert(ma);
  return xcalloc(nmemb, size);
}

void* ma_realloc(MA *ma, void *ptr, size_t size)
{
  assert(ma);
  return xrealloc(ptr, size);
}

void ma_free(void *ptr, MA *ma)
{
  assert(ma);
  free(ptr);
}

void ma_delete(MA *ma)
{
  if (!ma) return;
  free(ma);
}
