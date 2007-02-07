/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_STREAM_H
#define GENOME_STREAM_H

#include "genome.h"

typedef struct Genome_stream Genome_stream;

Genome_feature* genome_stream_next_feature(Genome_stream *gs);

genome_stream_push(Genome_stream *gs);
void genome_stream_free(Genome_stream *gs);

#endif
