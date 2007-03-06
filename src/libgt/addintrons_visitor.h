/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ADDINTRONS_VISITOR_H
#define ADDINTRONS_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct AddIntronsVisitor AddIntronsVisitor;

#include <libgt/genome_visitor.h>

const GenomeVisitorClass* addintrons_visitor_class(void);
GenomeVisitor*            addintrons_visitor_new(Env*);

#endif
