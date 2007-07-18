/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef UPGMA_H
#define UPGMA_H

#include <stdio.h>
#include "libgtcore/env.h"

typedef struct UPGMA UPGMA;

typedef double (*UPGMADistFunc)(unsigned long, unsigned long, void *data, Env*);

UPGMA* upgma_new(unsigned long num_of_taxa, void *data, UPGMADistFunc distfunc,
                 Env*);
void   upgma_show_tree(const UPGMA*, FILE*);
void   upgma_delete(UPGMA*, Env*);

#endif
