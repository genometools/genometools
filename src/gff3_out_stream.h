/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GFF3_OUT_STREAM_H
#define GFF3_OUT_STREAM_H

#include "genome_stream.h"

/* implements the ``genome stream'' interface */
typedef struct Gff3_out_stream Gff3_out_stream;

const Genome_stream_class* gff3_out_stream_class(void);
Genome_stream*             gff3_out_stream_new(Genome_stream*, FILE*);

#endif
