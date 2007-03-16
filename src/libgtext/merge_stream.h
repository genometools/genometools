/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MERGE_STREAM_H
#define MERGE_STREAM_H

#include <stdio.h>
#include <libgtext/genome_stream.h>

/* implements the ``genome_stream'' interface */
typedef struct MergeStream MergeStream;

const GenomeStreamClass* merge_stream_class(void);
GenomeStream*            merge_stream_new(const Array *genome_streams, Env*);

#endif
