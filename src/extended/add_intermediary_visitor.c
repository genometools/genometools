/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr.h"
#include "core/ma_api.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/add_intermediary_visitor.h"
#include "extended/node_visitor_rep.h"

struct GtAddIntermediaryVisitor {
  const GtNodeVisitor parent_instance;
  char *outside_type,
       *intermediary_type;
  GtFeatureNode *parent_feature,
                *previous_feature;
};

#define gt_add_intermediary_visitor_cast(GV)\
        gt_node_visitor_cast(gt_add_intermediary_visitor_class(), GV)

static int add_intermediary_in_children(GtGenomeNode *gn, void *data,
                                        GT_UNUSED GtError *err)
{
  GtAddIntermediaryVisitor *aiv = (GtAddIntermediaryVisitor*) data;
  GtFeatureNode *current_feature, *intermediary_node;
  GtRange previous_range, current_range, intermediary_range;
  GtStrand previous_strand, current_strand, intermediary_strand;
  GtStr *parent_seqid;
  gt_error_check(err);
  current_feature = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(current_feature);
  if (gt_feature_node_has_type(current_feature, aiv->outside_type)) {
    if (aiv->previous_feature) {
      /* determine intermediary range */
      previous_range = gt_genome_node_get_range((GtGenomeNode*)
                                                aiv->previous_feature);
      current_range = gt_genome_node_get_range(gn);
      gt_assert(previous_range.end < current_range.start);
      intermediary_range.start = previous_range.end + 1;
      intermediary_range.end = current_range.start - 1;

      /* determine intermediary strand */
      previous_strand = gt_feature_node_get_strand(aiv->previous_feature);
      current_strand = gt_feature_node_get_strand(current_feature);
      gt_assert(previous_strand == current_strand);
      intermediary_strand = previous_strand;

      /* determine sequence id */
      parent_seqid =
        gt_genome_node_get_seqid((GtGenomeNode*) aiv->parent_feature);
      gt_assert(!gt_str_cmp(parent_seqid,
                            gt_genome_node_get_seqid((GtGenomeNode*)
                                                     aiv->previous_feature)));
      gt_assert(!gt_str_cmp(parent_seqid,
                            gt_genome_node_get_seqid((GtGenomeNode*)
                                                     current_feature)));

      /* create intermediary */
      intermediary_node = (GtFeatureNode*)
                    gt_feature_node_new(parent_seqid, aiv->intermediary_type,
                                        intermediary_range.start,
                                        intermediary_range.end,
                                        intermediary_strand);
      gt_feature_node_add_child(aiv->parent_feature, intermediary_node);
    }
    aiv->previous_feature = current_feature;
  }
  return 0;
}

static int add_intermediary_if_necessary(GtGenomeNode *gn, void *data,
                                    GtError *err)
{
  GtAddIntermediaryVisitor *aiv = (GtAddIntermediaryVisitor*) data;
  GtFeatureNode *fn;
  gt_error_check(err);
  fn = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(fn);
  aiv->parent_feature = fn;
  aiv->previous_feature = NULL;
  return gt_genome_node_traverse_direct_children(gn, aiv,
                                                 add_intermediary_in_children,
                                                 err);
}

static void add_intermediary_visitor_free(GtNodeVisitor *nv)
{
  GtAddIntermediaryVisitor *aiv = gt_add_intermediary_visitor_cast(nv);
  gt_free(aiv->intermediary_type);
  gt_free(aiv->outside_type);
}

static int add_intermediary_visitor_feature_node(GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GtError *err)
{
  GtAddIntermediaryVisitor *v;
  gt_error_check(err);
  v = gt_add_intermediary_visitor_cast(nv);
  return gt_genome_node_traverse_children((GtGenomeNode*) fn, v,
                                          add_intermediary_if_necessary, false,
                                          err);
}

const GtNodeVisitorClass* gt_add_intermediary_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtAddIntermediaryVisitor),
                                    add_intermediary_visitor_free,
                                    NULL,
                                    add_intermediary_visitor_feature_node,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_add_intermediary_visitor_new(const char *outside_type,
                                               const char *intermediary_type)
{
  GtAddIntermediaryVisitor *aiv;
  GtNodeVisitor *nv;
  gt_assert(outside_type && intermediary_type);
  nv = gt_node_visitor_create(gt_add_intermediary_visitor_class());
  aiv = gt_add_intermediary_visitor_cast(nv);
  aiv->outside_type = gt_cstr_dup(outside_type);
  aiv->intermediary_type = gt_cstr_dup(intermediary_type);
  return nv;
}
