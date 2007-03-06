/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MUTATE_H
#define MUTATE_H

#include "alpha.h"
#include "seq.h"

/* returns a Seq which is the mutated sequence <orig_seq> of length <len> over
   alphabet <alpha>. <rate> denotes the error rate (must be >=0 && <= 100) */
Seq* mutate(const char *description, const char *orig_seq, unsigned long len,
            Alpha *alpha, unsigned int rate, Env*);

#endif
