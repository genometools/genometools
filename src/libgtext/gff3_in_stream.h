/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GFF3_IN_STREAM_H
#define GFF3_IN_STREAM_H

#include <stdio.h>
#include <libgtext/genome_stream.h>

/* implements the ``genome_stream'' interface */
typedef struct GFF3InStream GFF3InStream;

const GenomeStreamClass* gff3_in_stream_class(void);
void                     gff3_in_stream_set_offset(GenomeStream*, long);
GenomeStream*            gff3_in_stream_new_unsorted(int num_of_files,
                                                     const char **filenames,
                                                     bool be_verbose, Env*);
/* filename == NULL -> use stdin */
GenomeStream*            gff3_in_stream_new_sorted(const char *filename,
                                                   bool be_verbose, Env*);

#endif
