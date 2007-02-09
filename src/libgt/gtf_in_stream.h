/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTF_IN_STREAM_H
#define GTF_IN_STREAM_H

#include <stdio.h>
#include "genome_stream.h"

/* implements the ``genome_stream'' interface */
typedef struct Gtf_in_stream Gtf_in_stream;

const Genome_stream_class* gtf_in_stream_class(void);
/* filename == NULL -> use stdin */
Genome_stream*             gtf_in_stream_new(const char *filename,
                                             unsigned int be_tolerant);

#endif
