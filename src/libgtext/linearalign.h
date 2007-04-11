/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef LINEARALIGN_H
#define LINEARALIGN_H

#include <libgtext/alignment.h>

/* (globally) align <u> and <v> in linear space (unit cost) and return one
   optimal Alignment */
Alignment* linearalign(const char *u, unsigned long ulen,
                       const char *v, unsigned long vlen, Env*);

#endif
