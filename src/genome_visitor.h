/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_VISITOR_H
#define GENOME_VISITOR_H

#include "log.h"

/* the ``genome visitor'' interface, a visitor for genome nodes */
typedef struct Genome_visitor_class Genome_visitor_class;
typedef struct Genome_visitor Genome_visitor;

#include "comment.h"
#include "genome_feature.h"
#include "sequence_region.h"

void genome_visitor_visit_comment(Genome_visitor*, Comment*, Log*);
void genome_visitor_visit_genome_feature(Genome_visitor*, Genome_feature*,
                                         Log*);
void genome_visitor_visit_sequence_region(Genome_visitor*, Sequence_region*,
                                          Log*);
void genome_visitor_free(Genome_visitor *gv);

#endif
