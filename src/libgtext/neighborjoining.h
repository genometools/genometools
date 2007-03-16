/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef NEIGHBORJOINING_H
#define NEIGHBORJOINING_H

#include <stdio.h>
#include <gtcore.h>

typedef struct NeighborJoining NeighborJoining;

typedef double (*NeighborJoiningDistFunc)(unsigned long, unsigned long, void*,
                                          Env*);

NeighborJoining* neighborjoining_new(unsigned long num_of_taxa, void *data,
                                     NeighborJoiningDistFunc, Env*);
void             neighborjoining_show_tree(const NeighborJoining*, FILE*);
void             neighborjoining_delete(NeighborJoining*, Env*);

#endif
