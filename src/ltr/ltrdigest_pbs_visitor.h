/*
  Copyright (c) 2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef LTRDIGEST_PBS_VISITOR_H
#define LTRDIGEST_PBS_VISITOR_H

#include "core/bioseq.h"
#include "core/error_api.h"
#include "core/range_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping_api.h"

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtLTRdigestPBSVisitor GtLTRdigestPBSVisitor;

GtNodeVisitor* gt_ltrdigest_pbs_visitor_new(GtRegionMapping *rmap,
                                            unsigned int radius,
                                            unsigned int max_edist,
                                            GtRange alilen,
                                            GtRange offsetlen,
                                            GtRange trnaoffsetlen,
                                            int ali_score_match,
                                            int ali_score_mismatch,
                                            int ali_score_insertion,
                                            int ali_score_deletion,
                                            GtBioseq *trna_lib,
                                            GtError *err);

int            gt_ltrdigest_pbs_visitor_unit_test(GtError *err);

#endif
