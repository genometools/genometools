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
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/ma_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/feature_node.h"
#include "extended/inter_feature_visitor.h"
#include "extended/node_visitor_api.h"

struct GtInterFeatureVisitor {
  const GtNodeVisitor parent_instance;
  char *outside_type,
       *inter_type;
  GtFeatureNode *parent_feature,
                *previous_feature;
};

#define gt_inter_feature_visitor_cast(GV)\
        gt_node_visitor_cast(gt_inter_feature_visitor_class(), GV)

static int inter_feature_in_children(GtFeatureNode *current_feature, void *data,
                                     GT_UNUSED GtError *err)
{
  GtInterFeatureVisitor *aiv = (GtInterFeatureVisitor*) data;
  GtFeatureNode *inter_node;
  GtRange previous_range, current_range, inter_range;
  GtStrand previous_strand, /*current_strand, */inter_strand;
  GtStr *parent_seqid;
  gt_error_check(err);
  gt_assert(current_feature);
  if (gt_feature_node_has_type(current_feature, aiv->outside_type)) {
    if (aiv->previous_feature) {
      /* determine inter range */
      previous_range = gt_genome_node_get_range((GtGenomeNode*)
                                                aiv->previous_feature);
      current_range = gt_genome_node_get_range((GtGenomeNode*) current_feature);
      if (previous_range.end >= current_range.start) {
        gt_warning("overlapping boundary features %lu-%lu and %lu-%lu, "
                   "not placing '%s' inter-feature",
                   previous_range.start,
                   previous_range.end,
                   current_range.start,
                   current_range.end,
                   aiv->inter_type);
        return 0;
      }
      if (current_range.start - previous_range.end < 2) {
        gt_warning("no space for inter-feature '%s' between %lu and %lu",
                   aiv->inter_type,
                   previous_range.end,
                   current_range.start);
        return 0;
      }
      inter_range.start = previous_range.end + 1;
      inter_range.end = current_range.start - 1;

      /* determine inter strand */
      previous_strand = gt_feature_node_get_strand(aiv->previous_feature);
      /*current_strand = gt_feature_node_get_strand(current_feature);*/
      gt_assert(previous_strand == gt_feature_node_get_strand(current_feature));
      inter_strand = previous_strand;

      /* determine sequence id */
      parent_seqid =
        gt_genome_node_get_seqid((GtGenomeNode*) aiv->parent_feature);
      gt_assert(!gt_str_cmp(parent_seqid,
                            gt_genome_node_get_seqid((GtGenomeNode*)
                                                     aiv->previous_feature)));
      gt_assert(!gt_str_cmp(parent_seqid,
                            gt_genome_node_get_seqid((GtGenomeNode*)
                                                     current_feature)));

      /* create inter feature */
      inter_node = (GtFeatureNode*)
                   gt_feature_node_new(parent_seqid, aiv->inter_type,
                                       inter_range.start, inter_range.end,
                                       inter_strand);
      gt_feature_node_add_child(aiv->parent_feature, inter_node);
    }
    aiv->previous_feature = current_feature;
  }
  return 0;
}

static int inter_feature_if_necessary(GtFeatureNode *fn, void *data,
                                      GtError *err)
{
  GtInterFeatureVisitor *aiv = (GtInterFeatureVisitor*) data;
  gt_error_check(err);
  gt_assert(fn);
  aiv->parent_feature = fn;
  aiv->previous_feature = NULL;
  return gt_feature_node_traverse_direct_children(fn, aiv,
                                                  inter_feature_in_children,
                                                  err);
}

static void inter_feature_visitor_free(GtNodeVisitor *nv)
{
  GtInterFeatureVisitor *aiv = gt_inter_feature_visitor_cast(nv);
  gt_free(aiv->inter_type);
  gt_free(aiv->outside_type);
}

static int inter_feature_visitor_feature_node(GtNodeVisitor *nv,
                                              GtFeatureNode *fn, GtError *err)
{
  GtInterFeatureVisitor *v;
  gt_error_check(err);
  v = gt_inter_feature_visitor_cast(nv);
  return gt_feature_node_traverse_children(fn, v, inter_feature_if_necessary,
                                           false, err);
}

const GtNodeVisitorClass* gt_inter_feature_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtInterFeatureVisitor),
                                    inter_feature_visitor_free,
                                    NULL,
                                    inter_feature_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_inter_feature_visitor_new(const char *outside_type,
                                            const char *inter_type)
{
  GtInterFeatureVisitor *aiv;
  GtNodeVisitor *nv;
  gt_assert(outside_type && inter_type);
  nv = gt_node_visitor_create(gt_inter_feature_visitor_class());
  aiv = gt_inter_feature_visitor_cast(nv);
  aiv->outside_type = gt_cstr_dup(outside_type);
  aiv->inter_type = gt_cstr_dup(inter_type);
  return nv;
}
