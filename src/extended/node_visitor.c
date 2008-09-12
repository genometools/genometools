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

#include <stdlib.h>
#include "core/assert.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/node_visitor_rep.h"

GtNodeVisitor* gt_node_visitor_create(const GtNodeVisitorClass *gvc)
{
  GtNodeVisitor *gv;
  gt_assert(gvc && gvc->size);
  gv = gt_calloc(1, gvc->size);
  gv->c_class = gvc;
  return gv;
}

void* gt_node_visitor_cast(GT_UNUSED const GtNodeVisitorClass *gvc,
                          GtNodeVisitor *gv)
{
  gt_assert(gvc && gv && gv->c_class == gvc);
  return gv;
}

int gt_node_visitor_visit_comment_node(GtNodeVisitor *gv, GtCommentNode *cn,
                                      GtError *err)
{
  gt_error_check(err);
  gt_assert(gv && cn && gv->c_class);
  if (gv->c_class->comment_node)
    return gv->c_class->comment_node(gv, cn, err);
  return 0;
}

int gt_node_visitor_visit_feature_node(GtNodeVisitor *gv, GtFeatureNode *gf,
                                      GtError *err)
{
  gt_error_check(err);
  gt_assert(gv && gf && gv->c_class && gv->c_class->feature_node);
  return gv->c_class->feature_node(gv, gf, err);
}

int gt_node_visitor_visit_region_node(GtNodeVisitor *gv, GtRegionNode *rn,
                                     GtError *err)
{
  gt_error_check(err);
  gt_assert(gv && rn && gv->c_class);
  if (gv->c_class->region_node)
    return gv->c_class->region_node(gv, rn, err);
  return 0;
}

int gt_node_visitor_visit_sequence_node(GtNodeVisitor *gv, GtSequenceNode *sn,
                                       GtError *err)
{
  gt_error_check(err);
  gt_assert(gv && sn && gv->c_class);
  if (gv->c_class->sequence_node)
    return gv->c_class->sequence_node(gv, sn, err);
  return 0;
}

void gt_node_visitor_delete(GtNodeVisitor *gv)
{
  if (!gv) return;
  gt_assert(gv->c_class);
  if (gv->c_class->free)
    gv->c_class->free(gv);
  gt_free(gv);
}
