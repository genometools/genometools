/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

void warning(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);
  fprintf(stderr, "warning: ");
  (void) vfprintf(stderr, format, ap);
  (void) putc('\n', stderr);
  va_end(ap);
}
