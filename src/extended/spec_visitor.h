/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd.

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

#ifndef SPEC_VISITOR_H
#define SPEC_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtSpecVisitor GtSpecVisitor;

#include "core/error_api.h"
#include "extended/feature_index_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping_api.h"
#include "extended/spec_results.h"

GtNodeVisitor* gt_spec_visitor_new(const char *specfile, GtSpecResults *res,
                                   GtError *err);
void           gt_spec_visitor_report_runtime_errors(GtSpecVisitor *sv);
void           gt_spec_visitor_fail_on_runtime_error(GtSpecVisitor *sv);
void           gt_spec_visitor_add_feature_index(GtSpecVisitor *sv,
                                                 GtFeatureIndex *fi);
void           gt_spec_visitor_add_region_mapping(GtSpecVisitor *sv,
                                                  GtRegionMapping *rm);

#endif
