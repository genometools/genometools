/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_VISITOR_H
#define GENOME_VISITOR_H

/* the ``genome visitor'' interface, a visitor for genome nodes */
typedef struct GenomeVisitorClass GenomeVisitorClass;
typedef struct GenomeVisitor GenomeVisitor;

#include <libgtext/comment.h>
#include <libgtext/genome_feature.h>
#include <libgtext/sequence_region.h>

int  genome_visitor_visit_comment(GenomeVisitor*, Comment*, Env*);
int  genome_visitor_visit_genome_feature(GenomeVisitor*, GenomeFeature*, Env*);
int  genome_visitor_visit_sequence_region(GenomeVisitor*, SequenceRegion*,
                                          Env*);
void genome_visitor_delete(GenomeVisitor *gv, Env*);

#endif
