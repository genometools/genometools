/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef REGIONCOV_VISITOR_H
#define REGIONCOV_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct RegionCovVisitor RegionCovVisitor;

#include <libgtext/genome_visitor.h>

const GenomeVisitorClass* regioncov_visitor_class(void);
GenomeVisitor*            regioncov_visitor_new(unsigned long max_feature_dist,
                                                Env*);
void                      regioncov_visitor_show_coverage(GenomeVisitor*, Env*);

#endif
