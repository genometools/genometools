/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GFF3_VISITOR_H
#define GFF3_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct Gff3_visitor Gff3_visitor;

#include "genome_visitor.h"

const GenomeVisitorClass* gff3_visitor_class(void);
GenomeVisitor*             gff3_visitor_new(FILE*);

#endif
