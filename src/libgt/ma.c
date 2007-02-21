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

void ma_delete(MA *ma)
{
  if (!ma) return;
  free(ma);
}
