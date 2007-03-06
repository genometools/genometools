/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EXTRACTFEAT_VISITOR_H
#define EXTRACTFEAT_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct ExtractFeatVisitor ExtractFeatVisitor;

#include <stdbool.h>
#include <libgt/genome_visitor.h>
#include <libgt/regionmapping.h>

const GenomeVisitorClass* extractfeat_visitor_class(void);
GenomeVisitor*            extractfeat_visitor_new(RegionMapping *rm,
                                                  GenomeFeatureType, bool join,
                                                  bool translate, Env*);

#endif
