/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FILTER_VISITOR_H
#define FILTER_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct FilterVisitor FilterVisitor;

#include "genome_visitor.h"

const GenomeVisitorClass* filter_visitor_class(void);
GenomeVisitor*            filter_visitor_new(Str *seqid,
                                             unsigned long max_gene_length,
                                             double min_gene_score);
unsigned long             filter_visitor_node_buffer_size(GenomeVisitor*);
GenomeNode*               filter_visitor_get_node(GenomeVisitor*);

#endif
