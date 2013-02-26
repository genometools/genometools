/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/dlist.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_parser.h"
#include "extended/node_stream_api.h"
#include "extended/targetbest_select_stream.h"

struct GtTargetbestSelectStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtDlist *trees;
  GtDlistelem *next;
  GtHashmap *target_to_elem; /* maps the target ids to GtDlist elements */
  bool in_stream_processed;
};

#define targetbest_select_stream_cast(GS)\
        gt_node_stream_cast(gt_targetbest_select_stream_class(), GS);

static void build_key(GtStr *key, GtFeatureNode *feature, GtStr *target_id)
{
  gt_assert(key && feature && target_id);
  gt_str_reset(key);
  gt_str_append_str(key, gt_genome_node_get_seqid((GtGenomeNode*) feature));
  gt_str_append_char(key, '\t'); /* cannot occur in seqid or target_id */
  gt_str_append_str(key, target_id);
}

static void include_feature(GtDlist *trees, GtHashmap *target_to_elem,
                            GtFeatureNode *feature, GtStr *key)
{
  gt_dlist_add(trees, feature);
  gt_hashmap_add(target_to_elem, gt_cstr_dup(gt_str_get(key)),
              gt_dlist_last(trees));
}

static void remove_elem(GtDlistelem *elem, GtDlist *trees,
                        GtHashmap *target_to_elem, GtStr *key)
{
  GtGenomeNode *node = gt_dlistelem_get_data(elem);
  gt_genome_node_delete(node);
  gt_dlist_remove(trees, elem);
  gt_hashmap_remove(target_to_elem, gt_str_get(key));
}

static void replace_previous_elem(GtDlistelem *previous_elem,
                                  GtFeatureNode *current_feature,
                                  GtDlist *trees, GtHashmap *target_to_elem,
                                  GtStr *key)
{
  remove_elem(previous_elem, trees, target_to_elem, key);
  include_feature(trees, target_to_elem, current_feature, key);
}

static void select_targetbest(GtFeatureNode *current_feature,
                              GtDlist *trees, GtHashmap *target_to_elem)
{
  unsigned long num_of_targets;
  GtDlistelem *previous_elem;
  GtStr *first_target_id;
  const char *target;
  GT_UNUSED int had_err;
  gt_assert(current_feature && trees);
  target = gt_feature_node_get_attribute(current_feature, GT_GFF_TARGET);
  gt_assert(target);
  first_target_id = gt_str_new();
  had_err = gt_gff3_parser_parse_target_attributes(target, &num_of_targets,
                                                   first_target_id, NULL, NULL,
                                                   "", 0, NULL);
  gt_assert(!had_err);
  if (num_of_targets == 1) {
    GtStr *key = gt_str_new();
    build_key(key, current_feature, first_target_id);
    if (!(previous_elem = gt_hashmap_get(target_to_elem, gt_str_get(key)))) {
      /* element with this target_id not included yet -> include it */
      include_feature(trees, target_to_elem, current_feature, key);
    }
    else {
      GtFeatureNode *previous_feature = gt_dlistelem_get_data(previous_elem);
      /* element with this target_id included already -> compare them */
      if (gt_feature_node_get_score(current_feature) >
          gt_feature_node_get_score(previous_feature)) {
        /* current feature is better -> replace previous feature */
        replace_previous_elem(previous_elem, current_feature, trees,
                              target_to_elem, key);
      }
      else /* current feature is not better -> remove it */
        gt_genome_node_delete((GtGenomeNode*) current_feature);
    }
    gt_str_delete(key);
  }
  else
    gt_dlist_add(trees, current_feature);
  gt_str_delete(first_target_id);
}

static int targetbest_select_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                         GtError *err)
{
  GtTargetbestSelectStream *tfs;
  GtGenomeNode *node;
  int had_err = 0;
  gt_error_check(err);
  tfs = targetbest_select_stream_cast(ns);

  if (!tfs->in_stream_processed) {
    while (!(had_err = gt_node_stream_next(tfs->in_stream, &node, err)) &&
           node) {
      if (gt_feature_node_try_cast(node) &&
          gt_feature_node_get_attribute((GtFeatureNode*) node, "Target")) {
        select_targetbest((GtFeatureNode*) node, tfs->trees,
                          tfs->target_to_elem);
      }
      else
        gt_dlist_add(tfs->trees, node);
    }
    tfs->next = gt_dlist_first(tfs->trees);
    tfs->in_stream_processed = true;
  }

  if (!had_err) {
    gt_assert(tfs->in_stream_processed);
    if (tfs->next) {
      *gn = gt_dlistelem_get_data(tfs->next);
      tfs->next = gt_dlistelem_next(tfs->next);
    }
    else
      *gn = NULL;
    return 0;
  }

  return had_err;
}

static void targetbest_select_stream_free(GtNodeStream *ns)
{
  GtTargetbestSelectStream *tfs = targetbest_select_stream_cast(ns);
  for (; tfs->next != NULL; tfs->next = gt_dlistelem_next(tfs->next))
    gt_genome_node_delete(gt_dlistelem_get_data(tfs->next));
  gt_dlist_delete(tfs->trees);
  gt_hashmap_delete(tfs->target_to_elem);
  gt_node_stream_delete(tfs->in_stream);
}

const GtNodeStreamClass* gt_targetbest_select_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtTargetbestSelectStream),
                                   targetbest_select_stream_free,
                                   targetbest_select_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_targetbest_select_stream_new(GtNodeStream *in_stream)
{
  GtTargetbestSelectStream *tfs;
  GtNodeStream *ns;
  gt_assert(in_stream);
  ns = gt_node_stream_create(gt_targetbest_select_stream_class(),
                             gt_node_stream_is_sorted(in_stream));
  tfs = targetbest_select_stream_cast(ns);
  tfs->in_stream = gt_node_stream_ref(in_stream);
  tfs->in_stream_processed = false;
  tfs->trees = gt_dlist_new(NULL);
  tfs->target_to_elem = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
  return ns;
}
