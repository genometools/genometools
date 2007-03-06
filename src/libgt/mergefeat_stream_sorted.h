/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MERGEFEAT_STREAM_SORTED_H
#define MERGEFEAT_STREAM_SORTED_H

#include <stdio.h>
#include <libgt/genome_stream.h>

/* implements the ``genome_stream'' interface */
typedef struct MergefeatStreamSorted MergefeatStreamSorted;

const GenomeStreamClass* mergefeat_stream_sorted_class(void);
GenomeStream*            mergefeat_stream_sorted_new(GenomeStream*, Env*);

#endif
