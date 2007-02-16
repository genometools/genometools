/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_VISITOR_H
#define GENOME_VISITOR_H

#include "log.h"

/* the ``genome visitor'' interface, a visitor for genome nodes */
typedef struct GenomeVisitorClass GenomeVisitorClass;
typedef struct GenomeVisitor GenomeVisitor;

#include "comment.h"
#include "genome_feature.h"
#include "sequence_region.h"

int  genome_visitor_visit_comment(GenomeVisitor*, Comment*, Log*, Error*);
int  genome_visitor_visit_genome_feature(GenomeVisitor*, GenomeFeature*, Log*,
                                         Error*);
int  genome_visitor_visit_sequence_region(GenomeVisitor*, SequenceRegion*,
                                          Log*, Error*);
void genome_visitor_free(GenomeVisitor *gv);

#endif
