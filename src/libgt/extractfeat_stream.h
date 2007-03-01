/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EXTRACTFEAT_STREAM_H
#define EXTRACTFEAT_STREAM_H

#include <stdio.h>
#include "genome_stream.h"
#include "regionmapping.h"
#include "str.h"

/* implements the ``genome_stream'' interface */
typedef struct ExtractFeatStream ExtractFeatStream;

const GenomeStreamClass* extractfeat_stream_class(void);

/* create a ExtractFeatStream, takes ownership of RegionMapping  */
GenomeStream*            extractfeat_stream_new(GenomeStream*, RegionMapping*,
                                                GenomeFeatureType type,
                                                bool join, bool translate,
                                                Env*);

#endif
