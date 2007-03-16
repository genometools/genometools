/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef WARNING_H
#define WARNING_H

/* print a warning and continue */
void warning(const char *format, ...)
  __attribute__ ((format (printf, 1, 2)));

#endif
