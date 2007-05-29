/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CAIRO_STREAM_H
#define CAIRO_STREAM_H

#include <libgtext/genome_stream.h>

/* implements the ``genome stream'' interface */
typedef struct CairoStream CairoStream;

const GenomeStreamClass* cairo_stream_class(void);
GenomeStream*            cairo_stream_new(GenomeStream*, Str *seqid,
                                          unsigned long from, unsigned long to,
                                          const char *png_filename, int width,
                                          Env*);
void                     cairo_stream_draw(CairoStream*, bool verbose, Env*);

#endif
