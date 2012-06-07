/*
  Copyright (c) 2010-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/cstr_table_api.h"
#include "core/hashmap.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/queue_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/add_ids_visitor.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/gff3_defines.h"
#include "extended/node_visitor_api.h"

struct GtAddIDsVisitor {
  const GtNodeVisitor parent_instance;
  GtQueue *node_buffer;
  GtCstrTable *defined_seqids;
  GtHashmap *undefined_sequence_regions; /* contains all (automatically created)
                                            sequence regions */
  bool ensure_sorting;
};

#define add_ids_visitor_cast(GV)\
        gt_node_visitor_cast(gt_add_ids_visitor_class(), GV)

typedef struct {
  GtGenomeNode *sequence_region; /* the automatically created sequence region */
  GtArray *feature_nodes; /* the features nodes which belong to this region */
  bool is_circular;
} AutomaticSequenceRegion;

static AutomaticSequenceRegion* automatic_sequence_region_new(bool is_circular)
{
  AutomaticSequenceRegion *auto_sr;
  auto_sr = gt_malloc(sizeof (AutomaticSequenceRegion));
  auto_sr->feature_nodes = gt_array_new(sizeof (GtFeatureNode*));
  auto_sr->is_circular = is_circular;
  return auto_sr;
}

static void automatic_sequence_region_delete(AutomaticSequenceRegion *auto_sr)
{
  unsigned long i;
  if (!auto_sr) return;
  gt_genome_node_delete(auto_sr->sequence_region);
  for (i = 0; i < gt_array_size(auto_sr->feature_nodes); i++) {
    gt_genome_node_delete(*(GtGenomeNode**)
                              gt_array_get(auto_sr->feature_nodes, i));
  }
  gt_array_delete(auto_sr->feature_nodes);
  gt_free(auto_sr);
}

static void add_ids_visitor_free(GtNodeVisitor *nv)
{
  GtAddIDsVisitor *add_ids_visitor = add_ids_visitor_cast(nv);
  gt_hashmap_delete(add_ids_visitor->undefined_sequence_regions);
  gt_cstr_table_delete(add_ids_visitor->defined_seqids);
  gt_queue_delete(add_ids_visitor->node_buffer);
}

static int add_ids_visitor_comment_node(GtNodeVisitor *nv, GtCommentNode *c,
                                        GT_UNUSED GtError *err)
{
  GtAddIDsVisitor *add_ids_visitor;
  gt_error_check(err);
  add_ids_visitor = add_ids_visitor_cast(nv);
  gt_queue_add(add_ids_visitor->node_buffer, c);
  return 0;
}

static int add_ids_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                        GtError *err)
{
  AutomaticSequenceRegion *auto_sr;
  GtAddIDsVisitor *aiv;
  const char *seqid;
  bool is_circular;
  aiv = add_ids_visitor_cast(nv);
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) fn));
  if (aiv->ensure_sorting && !gt_cstr_table_get(aiv->defined_seqids, seqid)) {
    gt_error_set(err, "the file %s is not sorted (seqid \"%s\" on line %u has "
                 "not been previously introduced with a \"%s\" line)",
                 gt_genome_node_get_filename((GtGenomeNode*) fn), seqid,
                 gt_genome_node_get_line_number((GtGenomeNode*) fn),
                 GT_GFF_SEQUENCE_REGION);
    return -1;
  }
  if (!gt_cstr_table_get(aiv->defined_seqids, seqid)) {
    GtFeatureNodeIterator *fni;
    GtFeatureNode *node;
    GtRange range = gt_genome_node_get_range((GtGenomeNode*) fn);
    is_circular = gt_feature_node_get_attribute(fn, GT_GFF_IS_CIRCULAR)
                  ? true : false;
    if (!is_circular) {
      fni = gt_feature_node_iterator_new(fn);
      while ((node = gt_feature_node_iterator_next(fni))) {
        GtRange node_range = gt_genome_node_get_range((GtGenomeNode*) node);
        range = gt_range_join(&range, &node_range);
      }
      gt_feature_node_iterator_delete(fni);
    }
    /* sequence region has not been previously introduced -> check if one has
       already been created automatically */
    auto_sr = gt_hashmap_get(aiv->undefined_sequence_regions, seqid);
    if (!auto_sr) {
      GtStr *seqid_str;
      /* sequence region has not been createad automatically -> do it now */
      gt_warning("seqid \"%s\" on line %u in file \"%s\" has not been "
                 "previously introduced with a \"%s\" line, create such a line "
                 "automatically", seqid,
                 gt_genome_node_get_line_number((GtGenomeNode*) fn),
                 gt_genome_node_get_filename((GtGenomeNode*) fn),
                 GT_GFF_SEQUENCE_REGION);
      auto_sr = automatic_sequence_region_new(is_circular);
      seqid_str = gt_genome_node_get_seqid((GtGenomeNode*) fn);
      auto_sr->sequence_region = gt_region_node_new(seqid_str, range.start,
                                                               range.end);
      gt_hashmap_add(aiv->undefined_sequence_regions, gt_str_get(seqid_str),
                     auto_sr);
    }
    else {
      if (auto_sr->is_circular) {
        gt_assert(!is_circular); /* XXX */
      }
      else if (is_circular) {
        gt_assert(!auto_sr->is_circular); /* XXX */
        auto_sr->is_circular = true;
        gt_genome_node_set_range(auto_sr->sequence_region, &range);
      }
      else {
        GtRange joined_range,
                sr_range = gt_genome_node_get_range(auto_sr->sequence_region);
        /* update the range of the sequence region */
        joined_range = gt_range_join(&range, &sr_range);
        gt_genome_node_set_range(auto_sr->sequence_region, &joined_range);
      }
    }
    gt_array_add(auto_sr->feature_nodes, fn);
  }
  else
    gt_queue_add(aiv->node_buffer, fn);
  return 0;
}

static int add_ids_visitor_meta_node(GtNodeVisitor *nv, GtMetaNode *mn,
                                     GT_UNUSED GtError *err)
{
  GtAddIDsVisitor *add_ids_visitor;
  gt_error_check(err);
  add_ids_visitor = add_ids_visitor_cast(nv);
  gt_queue_add(add_ids_visitor->node_buffer, mn);
  return 0;
}

static int add_ids_visitor_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                       GT_UNUSED GtError *err)
{
  GtAddIDsVisitor *aiv;
  const char *seqid;
  int had_err = 0;
  gt_error_check(err);
  aiv = add_ids_visitor_cast(nv);
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) rn));
  if (gt_hashmap_get(aiv->undefined_sequence_regions, seqid)) {
    gt_error_set(err, "genome feature with id \"%s\" has been defined before "
                 "the corresponding \"%s\" definition on line %u in file "
                 "\"%s\"", seqid, GT_GFF_SEQUENCE_REGION,
                 gt_genome_node_get_line_number((GtGenomeNode*) rn),
                 gt_genome_node_get_filename((GtGenomeNode*) rn));
    had_err = -1;
  }
  if (!had_err) {
    if (!gt_cstr_table_get(aiv->defined_seqids, seqid))
      gt_cstr_table_add(aiv->defined_seqids, seqid);
    gt_queue_add(aiv->node_buffer, rn);
  }
  return had_err;
}

static int add_ids_visitor_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn,
                                         GT_UNUSED GtError *err)
{
  GtAddIDsVisitor *add_ids_visitor;
  gt_error_check(err);
  add_ids_visitor = add_ids_visitor_cast(nv);
  /* sequence nodes have to be at the end of a stream -> finalize first */
  gt_add_ids_visitor_finalize(nv);
  /* then add sequence node to buffer */
  gt_queue_add(add_ids_visitor->node_buffer, sn);
  return 0;
}

static int add_ids_visitor_eof_node(GtNodeVisitor *nv, GtEOFNode *eofn,
                                    GT_UNUSED GtError *err)
{
  GtAddIDsVisitor *add_ids_visitor;
  gt_error_check(err);
  add_ids_visitor = add_ids_visitor_cast(nv);
  gt_add_ids_visitor_finalize(nv);
  gt_queue_add(add_ids_visitor->node_buffer, eofn);
  return 0;
}

const GtNodeVisitorClass* gt_add_ids_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtAddIDsVisitor),
                                    add_ids_visitor_free,
                                    add_ids_visitor_comment_node,
                                    add_ids_visitor_feature_node,
                                    add_ids_visitor_region_node,
                                    add_ids_visitor_sequence_node,
                                    add_ids_visitor_eof_node);
    gt_node_visitor_class_set_meta_node_func(nvc, add_ids_visitor_meta_node);
  }
  return nvc;
}

GtNodeVisitor* gt_add_ids_visitor_new(bool ensure_sorting)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_add_ids_visitor_class());
  GtAddIDsVisitor *add_ids_visitor = add_ids_visitor_cast(nv);
  add_ids_visitor->node_buffer = gt_queue_new();
  add_ids_visitor->defined_seqids = gt_cstr_table_new();
  add_ids_visitor->undefined_sequence_regions =
    gt_hashmap_new(GT_HASH_STRING, NULL,
                   (GtFree) automatic_sequence_region_delete);
  add_ids_visitor->ensure_sorting = ensure_sorting;
  return nv;
}

unsigned long gt_add_ids_visitor_node_buffer_size(GtNodeVisitor *nv)
{
  GtAddIDsVisitor *add_ids_visitor = add_ids_visitor_cast(nv);
  return gt_queue_size(add_ids_visitor->node_buffer);
}

GtGenomeNode* gt_add_ids_visitor_get_node(GtNodeVisitor *nv)
{
  GtAddIDsVisitor *add_ids_visitor = add_ids_visitor_cast(nv);
  return gt_queue_get(add_ids_visitor->node_buffer);
}

static int add_auto_sr_to_queue(GT_UNUSED void *key, void *value, void *data,
                                GT_UNUSED GtError *err)
{
  AutomaticSequenceRegion *auto_sr = value;
  GtQueue *genome_nodes = data;
  GtGenomeNode *gf;
  unsigned int i;
  gt_error_check(err);
  gt_assert(key && value && data);
  if (gt_array_size(auto_sr->feature_nodes)) {
    gt_queue_add(genome_nodes, auto_sr->sequence_region);
    auto_sr->sequence_region = NULL;
    for (i = 0; i < gt_array_size(auto_sr->feature_nodes); i++) {
      gf = *(GtGenomeNode**) gt_array_get(auto_sr->feature_nodes, i);
      gt_queue_add(genome_nodes, gf);
    }
    gt_array_reset(auto_sr->feature_nodes);
  }
  return 0;
}

void gt_add_ids_visitor_finalize(GtNodeVisitor *nv)
{
  GT_UNUSED int had_err;
  GtAddIDsVisitor *add_ids_visitor = add_ids_visitor_cast(nv);
  had_err = gt_hashmap_foreach(add_ids_visitor->undefined_sequence_regions,
                               add_auto_sr_to_queue,
                               add_ids_visitor->node_buffer, NULL);
  gt_assert(!had_err); /* add_auto_sr_to_queue() is sane */
  gt_hashmap_reset(add_ids_visitor->undefined_sequence_regions);
}
