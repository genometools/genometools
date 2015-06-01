/*
  Copyright (c) 2012 Manuela Beckert <9beckert@informatik.uni-hamburg.de>
  Copyright (c) 2012 Dorle Osterode <9osterode@informatik.uni-hamburg.de>
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef TIR_STREAM_H
#define TIR_STREAM_H

#include "extended/node_stream_api.h"
#include "match/xdrop.h"

/* implements the "genome stream" interface */
typedef struct GtTIRStream GtTIRStream;

/* XXX: Needs documentation. */
const GtNodeStreamClass* gt_tir_stream_class(void);
GtNodeStream*            gt_tir_stream_new(GtStr *str_indexname,
                                           GtUword min_seed_length,
                                           GtUword min_TIR_length,
                                           GtUword max_TIR_length,
                                           GtUword min_TIR_distance,
                                           GtUword max_TIR_distance,
                                           GtXdropArbitraryscores arbit_scores,
                                           int xdrop_belowscore,
                                           double similarity_threshold,
                                           bool best_overlaps,
                                           bool longest_overlaps,
                                           bool no_overlaps,
                                           GtUword min_TSD_length,
                                           GtUword max_TSD_length,
                                           GtUword vicinity,
                                           GtError *err);

const GtEncseq* gt_tir_stream_get_encseq(GtTIRStream *ts);
void            gt_tir_stream_disable_md5_seqids(GtTIRStream *tir_stream);
void            gt_tir_stream_enable_md5_seqids(GtTIRStream *tir_stream);
void            gt_tir_stream_disable_seqids(GtTIRStream *tir_stream);
void            gt_tir_stream_enable_seqids(GtTIRStream *tir_stream);

#endif
