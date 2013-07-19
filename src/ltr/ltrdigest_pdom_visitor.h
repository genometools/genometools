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

#ifndef LTRDIGEST_PDOM_VISITOR_H
#define LTRDIGEST_PDOM_VISITOR_H

#include "extended/node_visitor.h"
#include "extended/region_mapping_api.h"
#include "ltr/pdom_model_set.h"

typedef enum {
  GT_PHMM_CUTOFF_TC,
  GT_PHMM_CUTOFF_GA,
  GT_PHMM_CUTOFF_NONE
} GtPdomCutoff;

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtLTRdigestPdomVisitor GtLTRdigestPdomVisitor;

GtNodeVisitor* gt_ltrdigest_pdom_visitor_new(GtPdomModelSet *model,
                                             double eval_cutoff,
                                             unsigned int chain_max_gap_length,
                                             GtPdomCutoff cutoff,
                                             GtRegionMapping *rmap,
                                             GtError *err);

void           gt_ltrdigest_pdom_visitor_output_all_chains(
                                                    GtLTRdigestPdomVisitor *lv);
void           gt_ltrdigest_pdom_visitor_set_root_type(
                                                     GtLTRdigestPdomVisitor *lv,
                                                     const char *type);
void           gt_ltrdigest_pdom_visitor_set_source_tag(
                                                     GtLTRdigestPdomVisitor *lv,
                                                     const char *tag);
#endif
