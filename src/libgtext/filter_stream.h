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

#ifndef FILTER_STREAM_H
#define FILTER_STREAM_H

#include "libgtext/genome_stream.h"

/* implements the ``genome stream'' interface */
typedef struct FilterStream FilterStream;

const GenomeStreamClass* filter_stream_class(void);
GenomeStream*            filter_stream_new(GenomeStream*,
                                           Str *seqid, Str *typefilter,
                                           Range contain_range,
                                           Range overlap_range, Strand strand,
                                           Strand targetstrand, bool has_CDS,
                                           unsigned long max_gene_length,
                                           unsigned long max_gene_num,
                                           double min_gene_score,
                                           double min_average_splice_site_prob);

#endif
