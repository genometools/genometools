/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GFF3_OUT_STREAM_H
#define GFF3_OUT_STREAM_H

#include <libgtext/genome_stream.h>

/* implements the ``genome stream'' interface */
typedef struct GFF3OutStream GFF3OutStream;

const GenomeStreamClass* gff3_out_stream_class(void);
GenomeStream*            gff3_out_stream_new(GenomeStream*, GenFile*, Env*);

#endif
