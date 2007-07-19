/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \file png_stream.h
 * \author Gordon Gremme <gremme@zbh.uni-hamburg.de>
 */

#ifndef PNG_STREAM_H
#define PNG_STREAM_H

#include "libgtext/genome_stream.h"

/* implements the ``genome stream'' interface */
typedef struct PNGStream PNGStream;

const GenomeStreamClass* png_stream_class(void);
GenomeStream*            png_stream_new(GenomeStream*, Str *seqid,
                                        unsigned long from, unsigned long to,
                                        const char *png_filename, int width,
                                        Env*);
void                     png_stream_draw(PNGStream*, bool verbose, Env*);

#endif
