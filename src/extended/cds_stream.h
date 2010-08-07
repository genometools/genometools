/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef CDS_STREAM_H
#define CDS_STREAM_H

#include <stdio.h>
#include "extended/node_stream_api.h"
#include "extended/region_mapping.h"

/* Implements the ``genome_stream'' interface. */
typedef struct GtCDSStream GtCDSStream;

const GtNodeStreamClass* gt_cds_stream_class(void);

/* Create a GtCDSStream, takes ownership of <region_mapping>. */
GtNodeStream*            gt_cds_stream_new(GtNodeStream *in_stream,
                                           GtRegionMapping *region_mapping,
                                           unsigned int minorflen,
                                           const char *source, bool start_codon,
                                           bool final_stop_codon,
                                           bool generic_start_codons);

#endif
