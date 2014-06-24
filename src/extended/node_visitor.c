/*
  Copyright (c) 2006-2012 Gordon Gremme <gordon@gremme.org>
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
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/node_visitor_api.h"

/* the <GtNodeVisitor> interface */
struct GtNodeVisitorClass {
  size_t size;
  GtNodeVisitorFreeFunc free;
  GtNodeVisitorCommentNodeFunc comment_node;
  GtNodeVisitorFeatureNodeFunc feature_node;
  GtNodeVisitorMetaNodeFunc meta_node;
  GtNodeVisitorRegionNodeFunc region_node;
  GtNodeVisitorSequenceNodeFunc sequence_node;
  GtNodeVisitorEOFNodeFunc eof_node;
};

GtNodeVisitorClass*
gt_node_visitor_class_new(size_t size,
                          GtNodeVisitorFreeFunc free,
                          GtNodeVisitorCommentNodeFunc comment_node,
                          GtNodeVisitorFeatureNodeFunc feature_node,
                          GtNodeVisitorRegionNodeFunc region_node,
                          GtNodeVisitorSequenceNodeFunc sequence_node,
                          GtNodeVisitorEOFNodeFunc eof_node)
{
  GtNodeVisitorClass *c_class;
  gt_assert(size);
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->comment_node = comment_node;
  c_class->feature_node = feature_node;
  c_class->region_node = region_node;
  c_class->sequence_node = sequence_node;
  c_class->eof_node = eof_node;
  return c_class;
}

void gt_node_visitor_class_set_meta_node_func(GtNodeVisitorClass *nvc,
                                              GtNodeVisitorMetaNodeFunc
                                              meta_node)
{
  gt_assert(nvc);
  nvc->meta_node = meta_node;
}

GtNodeVisitor* gt_node_visitor_create(const GtNodeVisitorClass *nvc)
{
  GtNodeVisitor *nv;
  gt_assert(nvc && nvc->size);
  nv = gt_calloc(1, nvc->size);
  nv->c_class = nvc;
  return nv;
}

void* gt_node_visitor_cast(GT_UNUSED const GtNodeVisitorClass *nvc,
                           GtNodeVisitor *nv)
{
  gt_assert(nvc && nv && nv->c_class == nvc);
  return nv;
}

int gt_node_visitor_visit_comment_node(GtNodeVisitor *nv, GtCommentNode *cn,
                                       GtError *err)
{
  gt_error_check(err);
  gt_assert(nv && cn && nv->c_class);
  if (nv->c_class->comment_node)
    return nv->c_class->comment_node(nv, cn, err);
  return 0;
}

int gt_node_visitor_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                       GtError *err)
{
  gt_error_check(err);
  gt_assert(nv && fn && nv->c_class);
  if (nv->c_class->feature_node)
    return nv->c_class->feature_node(nv, fn, err);
  return 0;
}

int gt_node_visitor_visit_meta_node(GtNodeVisitor *nv, GtMetaNode *mn,
                                    GtError *err)
{
  gt_error_check(err);
  gt_assert(nv && mn && nv->c_class);
  if (nv->c_class->meta_node)
    return nv->c_class->meta_node(nv, mn, err);
  return 0;
}

int gt_node_visitor_visit_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                      GtError *err)
{
  gt_error_check(err);
  gt_assert(nv && rn && nv->c_class);
  if (nv->c_class->region_node)
    return nv->c_class->region_node(nv, rn, err);
  return 0;
}

int gt_node_visitor_visit_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn,
                                        GtError *err)
{
  gt_error_check(err);
  gt_assert(nv && sn && nv->c_class);
  if (nv->c_class->sequence_node)
    return nv->c_class->sequence_node(nv, sn, err);
  return 0;
}

int gt_node_visitor_visit_eof_node(GtNodeVisitor *nv, GtEOFNode *eofn,
                                   GtError *err)
{
  gt_error_check(err);
  gt_assert(nv && eofn && nv->c_class);
  if (nv->c_class->eof_node)
    return nv->c_class->eof_node(nv, eofn, err);
  return 0;
}

void gt_node_visitor_delete(GtNodeVisitor *nv)
{
  if (!nv) return;
  gt_assert(nv->c_class);
  if (nv->c_class->free)
    nv->c_class->free(nv);
  gt_free(nv);
}
