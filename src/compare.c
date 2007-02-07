/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <string.h>
#include "compare.h"

int compare(const void *a, const void *b)
{
  return strcmp(*((const char**) a), *((const char**) b));
}
