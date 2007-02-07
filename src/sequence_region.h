/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SEQUENCE_REGION_H
#define SEQUENCE_REGION_H

/* implements the ``genome node'' interface */
typedef struct Sequence_region Sequence_region;

#include "genome_node.h"
#include "range.h"
#include "str.h"

const Genome_node_class* sequence_region_class(void);
Genome_node* sequence_region_new(Str *seqid,
                                 Range range,
                                 const char *filename,
                                 unsigned long line_number);

#endif
