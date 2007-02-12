/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FILTER_STREAM_H
#define FILTER_STREAM_H

#include "genome_stream.h"

/* implements the ``genome stream'' interface */
typedef struct Filter_stream Filter_stream;

const GenomeStreamClass* filter_stream_class(void);
GenomeStream*             filter_stream_new(GenomeStream*,
                                             unsigned long max_gene_length,
                                             double min_gene_score);

#endif
