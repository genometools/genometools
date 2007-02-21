/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef NEIGHBORJOINING_H
#define NEIGHBORJOINING_H

#include <stdio.h>

typedef struct NeighborJoining NeighborJoining;

NeighborJoining* neighborjoining_new(unsigned long num_of_taxa, void *data,
                                     double (*distfunc)
                                     (unsigned long, unsigned long, void*));
void             neighborjoining_show_tree(const NeighborJoining*, FILE*);
void             neighborjoining_delete(NeighborJoining*);

#endif
