/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ALIGN_H
#define ALIGN_H

#include "alignment.h"

/* (globally) align u and v (unit cost) and return one optimal Alignment */
Alignment* align(const char *u, unsigned long ulen,
                 const char *v, unsigned long vlen, Env*);

/* align u and v (unit cost), call proc_alignment for each optimal
   Alignment, and call proc_aligns with the number of optimal alignments*/
void align_all(const char *u, unsigned long ulen,
               const char *v, unsigned long vlen,
               void (*proc_alignment)(const Alignment*, void *data),
               void (*proc_aligns)(unsigned long, void *data), void *data,
               Env*);

#endif
