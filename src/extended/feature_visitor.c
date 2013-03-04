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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/unused_api.h"
#include "extended/feature_index.h"
#include "extended/feature_visitor.h"
#include "extended/node_visitor_api.h"
#include "extended/region_node.h"

struct GtFeatureVisitor {
  const GtNodeVisitor parent_instance;
        GtFeatureIndex *feature_index;
};

#define feature_visitor_cast(GV)\
        gt_node_visitor_cast(gt_feature_visitor_class(), GV)

static void feature_visitor_free(GtNodeVisitor *gv)
{
  GtFeatureVisitor *feature_visitor = feature_visitor_cast(gv);
  gt_assert(feature_visitor);
  gt_feature_index_delete(feature_visitor->feature_index);
}

static int feature_visitor_feature_node(GtNodeVisitor *gv, GtFeatureNode *fn,
                                        GtError *err)
{
  GtFeatureVisitor *v = feature_visitor_cast(gv);
  int rval = 0;
  gt_error_check(err);
  rval = gt_feature_index_add_feature_node(v->feature_index, fn, err);
  return rval;
}

static int feature_visitor_region_node(GtNodeVisitor *gv, GtRegionNode *rn,
                                       GtError *err)
{
  GtFeatureVisitor *v = feature_visitor_cast(gv);
  gt_error_check(err);
  return gt_feature_index_add_region_node(v->feature_index, rn, err);
}

const GtNodeVisitorClass* gt_feature_visitor_class()
{
  static const GtNodeVisitorClass *gvc = NULL;
  gt_class_alloc_lock_enter();
  if (!gvc) {
    gvc = gt_node_visitor_class_new(sizeof (GtFeatureVisitor),
                                    feature_visitor_free,
                                    NULL,
                                    feature_visitor_feature_node,
                                    feature_visitor_region_node,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return gvc;
}

GtNodeVisitor* gt_feature_visitor_new(GtFeatureIndex *fi)
{
  GtNodeVisitor *gv;
  GtFeatureVisitor *feature_visitor;
  gt_assert(fi != NULL);
  gv = gt_node_visitor_create(gt_feature_visitor_class());
  feature_visitor = feature_visitor_cast(gv);
  feature_visitor->feature_index = gt_feature_index_ref(fi);
  gt_assert(feature_visitor != NULL);
  return gv;
}
