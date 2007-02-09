/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MERGEFEAT_STREAM_SORTED_H
#define MERGEFEAT_STREAM_SORTED_H

#include <stdio.h>
#include "genome_stream.h"

/* implements the ``genome_stream'' interface */
typedef struct Mergefeat_stream_sorted Mergefeat_stream_sorted;

const Genome_stream_class* mergefeat_stream_sorted_class(void);
Genome_stream*             mergefeat_stream_sorted_new(Genome_stream*);

#endif
