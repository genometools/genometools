/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTF_VISITOR_H
#define GTF_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct GTFVisitor GTFVisitor;

#include <libgtext/genome_visitor.h>

const GenomeVisitorClass* gtf_visitor_class(void);
GenomeVisitor*            gtf_visitor_new(GenFile*, Env*);

#endif
