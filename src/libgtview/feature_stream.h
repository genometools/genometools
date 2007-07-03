/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FEATURE_STREAM_H
#define FEATURE_STREAM_H

#include <stdio.h>
#include <libgtext/genome_stream.h>
#include <libgtcore/hashtable.h>
#include <libgtview/feature_index.h>

/* implements the ``genome_stream'' interface */
typedef struct FeatureStream FeatureStream;

const GenomeStreamClass* feature_stream_class(void);

/* create a FeatureStream which writes to FeatureIndex */
GenomeStream*            feature_stream_new(GenomeStream*, FeatureIndex*, Env*);

#endif
