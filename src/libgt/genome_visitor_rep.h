/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_VISITOR_REP_H
#define GENOME_VISITOR_REP_H

#include <stdio.h>
#include "genome_visitor.h"

/* the ``genome visitor'' interface */
struct GenomeVisitorClass {
  size_t size;
  void (*free)(GenomeVisitor*, Env*);
  int  (*comment)(GenomeVisitor*, Comment*, Env*);
  int  (*genome_feature)(GenomeVisitor*, GenomeFeature*, Env*);
  int  (*sequence_region)(GenomeVisitor*, SequenceRegion*, Env*);
  int  (*default_func)(GenomeVisitor*, GenomeNode*, Env*);
};

struct GenomeVisitor {
  const GenomeVisitorClass *c_class;
};

GenomeVisitor* genome_visitor_create(const GenomeVisitorClass*, Env*);
void*          genome_visitor_cast(const GenomeVisitorClass*, GenomeVisitor*);

#endif
