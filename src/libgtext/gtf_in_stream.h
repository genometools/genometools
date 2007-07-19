/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTF_IN_STREAM_H
#define GTF_IN_STREAM_H

#include <stdio.h>
#include "libgtext/genome_stream.h"

/* implements the ``genome_stream'' interface */
typedef struct GTFInStream GTFInStream;

const GenomeStreamClass* gtf_in_stream_class(void);
/* filename == NULL -> use stdin */
GenomeStream*            gtf_in_stream_new(const char *filename,
                                           bool be_tolerant, Env*);

#endif
