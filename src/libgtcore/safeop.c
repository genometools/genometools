/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <libgtcore/safeop.h>

long safe_cast_to_long_type(unsigned long ulong, const char *filename,
                            unsigned int line)
{
  if (ulong > ~(1UL << ((8UL * sizeof (unsigned long)) - 1))) {
    fprintf(stderr, "%s, l.%u: illegal cast to long (ulong=%lu)\n",
            filename, line, ulong);
    exit(EXIT_FAILURE);
  }
  return (long) ulong;
}

unsigned long safe_cast_to_ulong_type(long slong, const char *filename,
                                      unsigned int line)
{
  if (slong < 0) {
    fprintf(stderr, "%s, l.%u: illegal cast to unsigned long (slong=%ld)\n",
            filename, line, slong);
    exit(EXIT_FAILURE);
  }
  return (unsigned long) slong;
}

