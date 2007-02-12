/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MERGEFEAT_STREAM_UNSORTED_H
#define MERGEFEAT_STREAM_UNSORTED_H

#include <stdio.h>
#include "genome_stream.h"

/* implements the ``genome_stream'' interface */
typedef struct Mergefeat_stream_unsorted Mergefeat_stream_unsorted;

const GenomeStreamClass* mergefeat_stream_unsorted_class(void);
GenomeStream*             mergefeat_stream_unsorted_new(GenomeStream*);

#endif
