/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef MD5_TO_SEQIDS_STREAM_H
#define MD5_TO_SEQIDS_STREAM_H

#include <stdio.h>
#include "extended/node_stream_api.h"
#include "extended/region_mapping.h"

/* Implements the ``genome_stream'' interface. */
typedef struct GtMD5ToSeqidsStream GtMD5ToSeqidsStream;

const GtNodeStreamClass* gt_md5_to_seqids_stream_class(void);

/* Create a GtSeqidToMD5Stream, takes ownership of <region_mapping>. */
GtNodeStream*            gt_md5_to_seqids_stream_new(GtNodeStream *in_stream,
                                                     GtRegionMapping
                                                     *region_mapping);

#endif
