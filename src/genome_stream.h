/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_STREAM_H
#define GENOME_STREAM_H

#include "genome_node.h"

/* the ``genome stream'' interface */
typedef struct Genome_stream_class Genome_stream_class;
typedef struct Genome_stream Genome_stream;

Genome_node* genome_stream_next_tree(Genome_stream*, Log*);
unsigned int genome_stream_is_sorted(Genome_stream*);
void         genome_stream_free(Genome_stream*);

#endif
