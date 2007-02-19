/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CDS_STREAM_H
#define CDS_STREAM_H

#include <stdio.h>
#include "genome_stream.h"

/* implements the ``genome_stream'' interface */
typedef struct CDS_stream CDS_stream;

const GenomeStreamClass* cds_stream_class(void);

GenomeStream*            cds_stream_new(GenomeStream*,
                                        const char *sequence_file,
                                        const char *source, Error*);

#endif
