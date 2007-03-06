/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GFF3_VISITOR_H
#define GFF3_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct GFF3Visitor GFF3Visitor;

#include <libgt/genome_visitor.h>

const GenomeVisitorClass* gff3_visitor_class(void);
GenomeVisitor*            gff3_visitor_new(FILE*, Env*);

#endif
