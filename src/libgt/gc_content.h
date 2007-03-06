/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GC_CONTENT_H
#define GC_CONTENT_H

#include <libgt/alpha.h>

/* show the GC-content for sequence <seq> with length <len> on stdout. <alpha>
   has to be compatible with a DNA alphabet */
void gc_content_show(const char *seq, unsigned long len, Alpha *alpha, Env*);

#endif
