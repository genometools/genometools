/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MERGEFEAT_VISITOR_H
#define MERGEFEAT_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct Mergefeat_visitor Mergefeat_visitor;

#include "genome_visitor.h"

const GenomeVisitorClass* mergefeat_visitor_class(void);
GenomeVisitor*             mergefeat_visitor_new(void);

#endif
