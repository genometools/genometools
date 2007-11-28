/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef STAT_STREAM_H
#define STAT_STREAM_H

#include <stdio.h>
#include "libgtext/genome_stream.h"

/* implements the ``genome_stream'' interface */
typedef struct StatStream StatStream;

const GenomeStreamClass* stat_stream_class(void);

GenomeStream*            stat_stream_new(GenomeStream*,
                                         bool gene_length_distri,
                                         bool gene_score_distri,
                                         bool exon_length_distri,
                                         bool exon_number_distri,
                                         bool intron_length_distri);
void                     stat_stream_show_stats(GenomeStream*);

#endif
