/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SORT_STREAM_H
#define SORT_STREAM_H

#include "genome_stream.h"

/* implements the ``genome stream'' interface */
typedef struct Sort_stream Sort_stream;

const Genome_stream_class* sort_stream_class(void);
Genome_stream*             sort_stream_new(Genome_stream*);

#endif
