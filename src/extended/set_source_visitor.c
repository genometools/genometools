/*
  Copyright (c) 2013 Daniel Standage <daniel.standage@gmail.com>

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.*/

#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/class_alloc_lock.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/feature_node_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/set_source_visitor.h"

#define set_source_visitor_cast(GV)\
        gt_node_visitor_cast(gt_set_source_visitor_class(), GV)

struct GtSetSourceVisitor {
  const GtNodeVisitor parent_instance;
  GtStr *newsource;
};

GtNodeVisitor *gt_set_source_visitor_new(GtStr *newsource)
{
  GtNodeVisitor *nv;
  gt_assert(newsource);
  nv = gt_node_visitor_create(gt_set_source_visitor_class());
  GtSetSourceVisitor *ssv = set_source_visitor_cast(nv);
  ssv->newsource = gt_str_ref(newsource);
  return nv;
}

static void set_source_visitor_free(GtNodeVisitor *nv)
{
  GtSetSourceVisitor *ssv = set_source_visitor_cast(nv);
  gt_str_delete(ssv->newsource);
}

static int set_source_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GT_UNUSED GtError *error)
{
  GtSetSourceVisitor *ssv;
  gt_error_check(error);
  ssv = set_source_visitor_cast(nv);

  GtFeatureNode *current;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  for (current  = gt_feature_node_iterator_next(iter);
       current != NULL;
       current  = gt_feature_node_iterator_next(iter))
  {
    gt_feature_node_set_source(current, ssv->newsource);
  }
  gt_feature_node_iterator_delete(iter);

  return 0;
}

const GtNodeVisitorClass *gt_set_source_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (GtSetSourceVisitor),
                                    set_source_visitor_free,
                                    NULL,
                                    set_source_visitor_visit_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}
