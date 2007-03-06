/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CDS_STREAM_H
#define CDS_STREAM_H

#include <stdio.h>
#include <libgt/genome_stream.h>
#include <libgt/regionmapping.h>

/* implements the ``genome_stream'' interface */
typedef struct CDSStream CDSStream;

const GenomeStreamClass* cds_stream_class(void);

/* create a CDSSTream, takes ownership of RegionMapping */
GenomeStream*            cds_stream_new(GenomeStream*, RegionMapping*,
                                        const char *source, Env*);

#endif
