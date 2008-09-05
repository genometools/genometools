/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef GFF3_IN_STREAM_H
#define GFF3_IN_STREAM_H

#include <stdio.h>
#include "extended/genome_stream.h"

/* implements the ``genome_stream'' interface */
typedef struct GFF3InStream GFF3InStream;

const GenomeStreamClass* gff3_in_stream_class(void);
/* Does not take ownership of <feature_type_factory>. */
void                     gff3_in_stream_set_feature_type_factory(GenomeStream*,
                                                         FeatureTypeFactory
                                                         *feature_type_factory);
/* Returns a StrArray which contains all type names in alphabetical order which
   have been parsed by this stream. The caller is responsible to free it! */
StrArray*                gff3_in_stream_get_used_types(GenomeStream*);
void                     gff3_in_stream_set_offset(GenomeStream*, long);
int                      gff3_in_stream_set_offsetfile(GenomeStream*, Str*,
                                                       Error*);
void                     gff3_in_stream_enable_tidy_mode(GenomeStream*);
GenomeStream*            gff3_in_stream_new_unsorted(int num_of_files,
                                                     const char **filenames,
                                                     bool be_verbose,
                                                     bool checkids);
/* filename == NULL -> use stdin */
GenomeStream*            gff3_in_stream_new_sorted(const char *filename,
                                                   bool be_verbose);

#endif
