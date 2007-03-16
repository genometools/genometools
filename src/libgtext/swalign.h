/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SWALIGN_H
#define SWALIGN_H

#include <libgtext/alignment.h>

/* (locally) align u and v (Smith-Waterman algorithm ) with the given score
   function and return one optimal Alignment.
   If no such alignment was found, NULL is returned. */
Alignment* swalign(Seq *u, Seq *v, const ScoreFunction*, Env*);

#endif
