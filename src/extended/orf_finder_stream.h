/*
  Copyright (c) 2010      Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef ORF_FINDER_STREAM_H
#define ORF_FINDER_STREAM_H

#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/hashmap_api.h"
#include "extended/node_stream_api.h"
#include "extended/orf_iterator_api.h"

/* implements the ``node stream'' interface */
typedef struct GtORFFinderStream GtORFFinderStream;

const GtNodeStreamClass* gt_orf_finder_stream_class(void);

GtNodeStream* gt_orf_finder_stream_new(GtNodeStream *in_stream,
                                       GtEncseq *encseq,
                                       GtHashmap *types,
                                       unsigned int min,
                                       unsigned int max,
                                       bool all,
                                       GtError *err);

#endif
