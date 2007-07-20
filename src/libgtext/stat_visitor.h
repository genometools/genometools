/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STAT_VISITOR_H
#define STAT_VISITOR_H

/* implements the ``genome visitor'' interface, gathers statistics */
typedef struct StatVisitor StatVisitor;

#include "libgtext/genome_visitor.h"

const GenomeVisitorClass* stat_visitor_class(void);
GenomeVisitor*            stat_visitor_new(bool gene_length_distri,
                                           bool gene_score_distri,
                                           bool exon_length_distri,
                                           bool intron_length_distri, Env*);
void                      stat_visitor_show_stats(GenomeVisitor*);

#endif
