/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTF_OUT_STREAM_H
#define GTF_OUT_STREAM_H

#include <libgtext/genome_stream.h>

/* implements the ``genome stream'' interface */
typedef struct GTFOutStream GTFOutStream;

const GenomeStreamClass* gtf_out_stream_class(void);
GenomeStream*            gtf_out_stream_new(GenomeStream*, GenFile*, Env*);

#endif
