/*
  Copyright (c) 2010-2011 Gordon Gremme <gordon@gremme.org>

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

#ifndef BSSM_TRAIN_VISITOR_H
#define BSSM_TRAIN_VISITOR_H

/* Implements the ``genome visitor'' interface. */
typedef struct GthBSSMTrainVisitor GthBSSMTrainVisitor;

#include "extended/node_visitor.h"
#include "extended/region_mapping_api.h"
#include "gth/bssm_seq_processor.h"

const GtNodeVisitorClass* gth_bssm_train_visitor_class(void);
GtNodeVisitor*            gth_bssm_train_visitor_new(GtRegionMapping
                                                     *region_mapping,
                                                     GthBSSMSeqProcessor*,
                                                     const char *filter_type,
                                                     const char *extract_type,
                                                     unsigned int
                                                     good_exon_count,
                                                     double cutoff);

#endif
