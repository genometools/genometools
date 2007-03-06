/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FILTER_STREAM_H
#define FILTER_STREAM_H

#include <libgt/genome_stream.h>

/* implements the ``genome stream'' interface */
typedef struct FilterStream FilterStream;

const GenomeStreamClass* filter_stream_class(void);
GenomeStream*            filter_stream_new(GenomeStream*,
                                           Str *seqid, Str *typefilter,
                                           unsigned long max_gene_length,
                                           unsigned long max_gene_num,
                                           double min_gene_score, Env*);

#endif
