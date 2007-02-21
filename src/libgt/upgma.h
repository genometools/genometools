/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef UPGMA_H
#define UPGMA_H

#include <stdio.h>

typedef struct UPGMA UPGMA;

UPGMA* upgma_new(unsigned long num_of_taxa, void *data,
                 double (*distfunc)(unsigned long, unsigned long, void *data));
void   upgma_show_tree(const UPGMA*, FILE*);
void   upgma_delete(UPGMA*);

#endif
