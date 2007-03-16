/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef AFFINEALIGN_H
#define AFFINEALIGN_H

#include <libgtext/alignment.h>

/* (globally) align u and v (affine gap costs) and return one optimal
   Alignment */
Alignment* affinealign(const char *u, unsigned long ulen,
                       const char *v, unsigned long vlen, int replacement_cost,
                       int gap_opening_cost, int gap_extension_cost, Env*);

#endif
