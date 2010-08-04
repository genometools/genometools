/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/undef.h"
#include "extended/seqids_to_md5_visitor.h"
#include "extended/node_visitor_rep.h"

struct GtSeqidsToMD5Visitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *region_mapping;
};

#define  seqids_to_md5_visitor_cast(GV)\
         gt_node_visitor_cast(gt_seqids_to_md5_visitor_class(), GV)

static void seqids_to_md5_visitor_free(GtNodeVisitor *nv)
{
  GtSeqidsToMD5Visitor *seqids_to_md5_visitor = seqids_to_md5_visitor_cast(nv);
  gt_assert(seqids_to_md5_visitor);
  gt_region_mapping_delete(seqids_to_md5_visitor->region_mapping);
}

static int seqids_to_md5_visitor_feature_node(GtNodeVisitor *nv,
                                              GtFeatureNode *fn,
                                              GtError *err)
{
  GtSeqidsToMD5Visitor *v = seqids_to_md5_visitor_cast(nv);
  gt_error_check(err);

  /* XXX */
  gt_assert(v && fn);

  return 0;
}

static int seqids_to_md5_visitor_region_node(GtNodeVisitor *nv,
                                             GtRegionNode *rn,
                                             GtError *err)
{
  GtSeqidsToMD5Visitor *v = seqids_to_md5_visitor_cast(nv);
  gt_error_check(err);

  /* XXX */
  gt_assert(v && rn);

  return 0;
}

const GtNodeVisitorClass* gt_seqids_to_md5_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtSeqidsToMD5Visitor),
                                    seqids_to_md5_visitor_free,
                                    NULL,
                                    seqids_to_md5_visitor_feature_node,
                                    seqids_to_md5_visitor_region_node,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_seqids_to_md5_visitor_new(GtRegionMapping *region_mapping)
{
  GtNodeVisitor *nv;
  GtSeqidsToMD5Visitor *seqids_to_md5_visitor;
  nv = gt_node_visitor_create(gt_seqids_to_md5_visitor_class());
  seqids_to_md5_visitor = seqids_to_md5_visitor_cast(nv);
  seqids_to_md5_visitor->region_mapping = region_mapping;
  return nv;
}
