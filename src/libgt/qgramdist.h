/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef QGRAMDIST_H
#define QGRAMDIST_H

#include "seq.h"

/* returns the q-gram distance of two Seqs. The alphabets of the given Seqs have
   to be compatible */
unsigned long qgramdist(Seq*, Seq*, unsigned int q, Env*);

#endif
