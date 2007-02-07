/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_VISITOR_REP_H
#define GENOME_VISITOR_REP_H

#include <stdio.h>
#include "genome_visitor.h"

/* the ``genome visitor'' interface */
struct Genome_visitor_class {
  size_t size;
  void (*free)(Genome_visitor*);
  void (*comment)(Genome_visitor*, Comment*, Log*);
  void (*genome_feature)(Genome_visitor*, Genome_feature*, Log*);
  void (*sequence_region)(Genome_visitor*, Sequence_region*, Log*);
  void (*default_func)(Genome_visitor*, Genome_node*, Log*);
};

struct Genome_visitor {
  const Genome_visitor_class *c_class;
};

Genome_visitor* genome_visitor_create(const Genome_visitor_class*);
void*           genome_visitor_cast(const Genome_visitor_class*,
                                    Genome_visitor*);

#endif
