/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CDS_VISITOR_H
#define CDS_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct CDSVisitor CDSVisitor;

#include "genome_visitor.h"

const GenomeVisitorClass* cds_visitor_class(void);
GenomeVisitor*            cds_visitor_new(Str *sequence_file, Str *source);

#endif
