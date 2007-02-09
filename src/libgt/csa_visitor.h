/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CSA_VISITOR_H
#define CSA_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct Csa_visitor Csa_visitor;

#include "genome_visitor.h"

const Genome_visitor_class* csa_visitor_class(void);
Genome_visitor*             csa_visitor_new(unsigned long join_length);
unsigned long               csa_visitor_node_buffer_size(Genome_visitor*);
Genome_node*                csa_visitor_get_node(Genome_visitor*);
void                        csa_visitor_process_cluster(Genome_visitor*,
                                                        unsigned int
                                                        final_cluster, Log*);

#endif
