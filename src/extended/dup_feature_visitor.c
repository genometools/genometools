/*
  Copyright (c) 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "extended/dup_feature_visitor.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"

struct GtDupFeatureVisitor {
  const GtNodeVisitor parent_instance;
  char *dest_type,
       *source_type;
};

#define gt_dup_feature_visitor_cast(GV)\
        gt_node_visitor_cast(gt_dup_feature_visitor_class(), GV)

static void dup_feature_visitor_free(GtNodeVisitor *nv)
{
  GtDupFeatureVisitor *aiv = gt_dup_feature_visitor_cast(nv);
  gt_free(aiv->source_type);
  gt_free(aiv->dest_type);
}

static GtFeatureNode* duplicate_feature(const GtFeatureNode *fn,
                                        const char *dest_type)
{
  GtFeatureNode *dup;
  GtStrArray *attr;
  unsigned long i;
  dup = (GtFeatureNode*)
        gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*) fn),
                            dest_type,
                            gt_genome_node_get_start((GtGenomeNode*) fn),
                            gt_genome_node_get_end((GtGenomeNode*) fn),
                            gt_feature_node_get_strand(fn));
  if (gt_feature_node_has_source(fn)) {
    GtStr *source = gt_str_new_cstr(gt_feature_node_get_source(fn));
    gt_feature_node_set_source(dup, source);
    gt_str_delete(source);
  }
  if (gt_feature_node_score_is_defined(fn))
    gt_feature_node_set_score(dup, gt_feature_node_get_score(fn));
  gt_feature_node_set_phase(dup, gt_feature_node_get_phase(fn));
  attr = gt_feature_node_get_attribute_list(fn);
  for (i = 0; i < gt_str_array_size(attr); i++) {
    const char *tag = gt_str_array_get(attr, i);
    gt_feature_node_add_attribute(dup, tag,
                                  gt_feature_node_get_attribute(fn, tag));
  }
  gt_str_array_delete(attr);
  return dup;
}

static int dup_feature_visitor_feature_node(GtNodeVisitor *nv,
                                            GtFeatureNode *fn,
                                            GT_UNUSED GtError *err)
{
  GtFeatureNodeIterator *fni, *child_iterator;
  GtFeatureNode *node, *child;
  GtDupFeatureVisitor *v;
  gt_error_check(err);
  v = gt_dup_feature_visitor_cast(nv);
  fni = gt_feature_node_iterator_new(fn);
  while ((node = gt_feature_node_iterator_next(fni))) {
    if (node != fn) {
      child_iterator = gt_feature_node_iterator_new_direct(node);
      while ((child = gt_feature_node_iterator_next(child_iterator))) {
        if (gt_feature_node_has_type(child, v->source_type)) {
          gt_feature_node_add_child(node,
                                    duplicate_feature(child,v->dest_type));
        }
      }
      gt_feature_node_iterator_delete(child_iterator);
    }
  }
  gt_feature_node_iterator_delete(fni);
  return 0;
}

const GtNodeVisitorClass* gt_dup_feature_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtDupFeatureVisitor),
                                    dup_feature_visitor_free,
                                    NULL,
                                    dup_feature_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_dup_feature_visitor_new(const char *dest_type,
                                          const char *source_type)
{
  GtDupFeatureVisitor *aiv;
  GtNodeVisitor *nv;
  gt_assert(dest_type && source_type);
  nv = gt_node_visitor_create(gt_dup_feature_visitor_class());
  aiv = gt_dup_feature_visitor_cast(nv);
  aiv->dest_type = gt_cstr_dup(dest_type);
  aiv->source_type = gt_cstr_dup(source_type);
  return nv;
}
