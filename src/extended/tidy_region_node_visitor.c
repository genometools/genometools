/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/queue.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/tidy_region_node_visitor.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"

struct GtTidyRegionNodeVisitor {
  const GtNodeVisitor parent_instance;
  GtQueue *node_buffer;
  GtHashmap *region_nodes;
};

#define tidy_region_node_visitor_cast(GV)\
        gt_node_visitor_cast(gt_tidy_region_node_visitor_class(), GV)

static void tidy_region_node_visitor_free(GtNodeVisitor *nv)
{
  GtTidyRegionNodeVisitor *tidy_region_node_visitor =
                                              tidy_region_node_visitor_cast(nv);
  gt_hashmap_delete(tidy_region_node_visitor->region_nodes);
  gt_queue_delete(tidy_region_node_visitor->node_buffer);
}

static int tidy_region_node_visitor_feature_node(GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GT_UNUSED GtError *err)
{
  GtTidyRegionNodeVisitor *trnv;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  GtGenomeNode *sr;
  GtRange range = {0, 0},
          sr_range,
          joined_range = {0, 0};
  bool first = true;
  int had_err = 0;
  const char *seqid;
  trnv = tidy_region_node_visitor_cast(nv);
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) fn));

  fni = gt_feature_node_iterator_new(fn);
  while ((node = gt_feature_node_iterator_next(fni))) {
    GtRange node_range = gt_genome_node_get_range((GtGenomeNode*) node);
    if (first) {
      range = node_range;
      first = false;
    }
    else
      range = gt_range_join(&range, &node_range);
  }
  gt_feature_node_iterator_delete(fni);

  sr = gt_hashmap_get(trnv->region_nodes, seqid);
  if (!sr) {
    gt_error_set(err, "seqid '%s' on line %u in file \"%s\" has not been "
                      "defined yet",
                 seqid,
                 gt_genome_node_get_line_number((GtGenomeNode*) fn),
                 gt_genome_node_get_filename((GtGenomeNode*) fn));
    had_err = -1;
  }

  if (!had_err) {
    gt_assert(sr);
    sr_range = gt_genome_node_get_range(sr);
    joined_range = gt_range_join(&range, &sr_range);
    gt_genome_node_set_range(sr, &joined_range);
    gt_queue_add(trnv->node_buffer, fn);
  }

  return 0;
}

static int tidy_region_node_visitor_region_node(GtNodeVisitor *nv,
                                                GtRegionNode *rn,
                                                GT_UNUSED GtError *err)
{
  GtTidyRegionNodeVisitor *trnv;
  const char *seqid;
  gt_error_check(err);
  trnv = tidy_region_node_visitor_cast(nv);
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) rn));
  gt_assert(seqid);
  if (!gt_hashmap_get(trnv->region_nodes, seqid)) {
    gt_hashmap_add(trnv->region_nodes, gt_cstr_dup(seqid),
                   gt_genome_node_ref((GtGenomeNode*) rn));
    gt_assert(gt_hashmap_get(trnv->region_nodes, seqid));
  }
  gt_queue_add(trnv->node_buffer, rn);
  return 0;
}

static int tidy_region_node_visitor_sequence_node(GtNodeVisitor *nv,
                                                  GtSequenceNode *sn,
                                                  GT_UNUSED GtError *err)
{
  GtTidyRegionNodeVisitor *trnv;
  gt_error_check(err);
  trnv = tidy_region_node_visitor_cast(nv);
  gt_queue_add(trnv->node_buffer, sn);
  return 0;
}

static int tidy_region_node_visitor_comment_node(GtNodeVisitor *nv,
                                                 GtCommentNode *cn,
                                                 GT_UNUSED GtError *err)
{
  GtTidyRegionNodeVisitor *trnv;
  gt_error_check(err);
  trnv = tidy_region_node_visitor_cast(nv);
  gt_queue_add(trnv->node_buffer, cn);
  return 0;
}

static int tidy_region_node_visitor_meta_node(GtNodeVisitor *nv,
                                              GtMetaNode *mn,
                                              GT_UNUSED GtError *err)
{
  GtTidyRegionNodeVisitor *trnv;
  gt_error_check(err);
  trnv = tidy_region_node_visitor_cast(nv);
  gt_queue_add(trnv->node_buffer, mn);
  return 0;
}

static int tidy_region_node_visitor_eof_node(GtNodeVisitor *nv,
                                             GtEOFNode *en,
                                             GT_UNUSED GtError *err)
{
  GtTidyRegionNodeVisitor *trnv;
  gt_error_check(err);
  trnv = tidy_region_node_visitor_cast(nv);
  gt_queue_add(trnv->node_buffer, en);
  return 0;
}

const GtNodeVisitorClass* gt_tidy_region_node_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtTidyRegionNodeVisitor),
                                    tidy_region_node_visitor_free,
                                    tidy_region_node_visitor_comment_node,
                                    tidy_region_node_visitor_feature_node,
                                    tidy_region_node_visitor_region_node,
                                    tidy_region_node_visitor_sequence_node,
                                    tidy_region_node_visitor_eof_node);
    gt_node_visitor_class_set_meta_node_func(nvc,
                                            tidy_region_node_visitor_meta_node);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_tidy_region_node_visitor_new(void)
{
  GtNodeVisitor *nv;
  GtTidyRegionNodeVisitor *tidy_region_node_visitor;
  nv = gt_node_visitor_create(gt_tidy_region_node_visitor_class());
  tidy_region_node_visitor = tidy_region_node_visitor_cast(nv);
  tidy_region_node_visitor->node_buffer = gt_queue_new();
  tidy_region_node_visitor->region_nodes = gt_hashmap_new(GT_HASH_STRING,
                                                gt_free_func,
                                                (GtFree) gt_genome_node_delete);
  return nv;
}

unsigned long gt_tidy_region_node_visitor_node_buffer_size(GtNodeVisitor *nv)
{
  GtTidyRegionNodeVisitor *tidy_region_node_visitor =
                                              tidy_region_node_visitor_cast(nv);
  return gt_queue_size(tidy_region_node_visitor->node_buffer);
}

GtGenomeNode* gt_tidy_region_node_visitor_get_node(GtNodeVisitor *nv)
{
  GtTidyRegionNodeVisitor *tidy_region_node_visitor =
                                              tidy_region_node_visitor_cast(nv);
  return gt_queue_get(tidy_region_node_visitor->node_buffer);
}
