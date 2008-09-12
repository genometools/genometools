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

#include <assert.h>
#include "core/cstr.h"
#include "core/dlist.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/undef.h"
#include "extended/genome_stream_rep.h"
#include "extended/gff3_parser.h"
#include "extended/targetbest_filter_stream.h"

struct TargetbestFilterStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GtDlist *trees;
  GtDlistelem *next;
  Hashmap *target_to_elem; /* maps the target ids to GtDlist elements */
  bool in_stream_processed;
};

#define targetbest_filter_stream_cast(GS)\
        genome_stream_cast(targetbest_filter_stream_class(), GS);

static void build_key(GtStr *key, GtFeatureNode *feature, GtStr *target_id)
{
  assert(key && feature && target_id);
  gt_str_reset(key);
  gt_str_append_str(key, gt_genome_node_get_seqid((GtGenomeNode*) feature));
  gt_str_append_char(key, '\t'); /* cannot occur in seqid or target_id */
  gt_str_append_str(key, target_id);
}

static void include_feature(GtDlist *trees, Hashmap *target_to_elem,
                            GtFeatureNode *feature, GtStr *key)
{
  gt_dlist_add(trees, feature);
  hashmap_add(target_to_elem, gt_cstr_dup(gt_str_get(key)),
              gt_dlist_last(trees));
}

static void remove_elem(GtDlistelem *elem, GtDlist *trees,
                        Hashmap *target_to_elem, GtStr *key)
{
  GtGenomeNode *node = gt_dlistelem_get_data(elem);
  gt_genome_node_rec_delete(node);
  gt_dlist_remove(trees, elem);
  hashmap_remove(target_to_elem, gt_str_get(key));
}

static void replace_previous_elem(GtDlistelem *previous_elem,
                                  GtFeatureNode *current_feature,
                                  GtDlist *trees, Hashmap *target_to_elem,
                                  GtStr *key)
{
  remove_elem(previous_elem, trees, target_to_elem, key);
  include_feature(trees, target_to_elem, current_feature, key);
}

static void filter_targetbest(GtFeatureNode *current_feature,
                              GtDlist *trees, Hashmap *target_to_elem)
{
  unsigned long num_of_targets;
  GtDlistelem *previous_elem;
  GtStr *first_target_id;
  const char *target;
  int had_err;
  assert(current_feature && trees);
  target = gt_feature_node_get_attribute((GtGenomeNode*) current_feature,
                                        TARGET_STRING);
  assert(target);
  first_target_id = gt_str_new();
  had_err = gt_gff3_parser_parse_target_attributes(target, &num_of_targets,
                                                   first_target_id, NULL, NULL,
                                                   "", 0, NULL);
  assert(!had_err);
  if (num_of_targets == 1) {
    GtStr *key = gt_str_new();
    build_key(key, current_feature, first_target_id);
    if (!(previous_elem = hashmap_get(target_to_elem, gt_str_get(key)))) {
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
        gt_genome_node_rec_delete((GtGenomeNode*) current_feature);
    }
    gt_str_delete(key);
  }
  else
    gt_dlist_add(trees, current_feature);
  gt_str_delete(first_target_id);
}

static int targetbest_filter_stream_next_tree(GenomeStream *gs,
                                              GtGenomeNode **gn, GtError *err)
{
  TargetbestFilterStream *tfs;
  GtGenomeNode *node;
  int had_err = 0;
  gt_error_check(err);
  tfs = targetbest_filter_stream_cast(gs);

  if (!tfs->in_stream_processed) {
    while (!(had_err = genome_stream_next_tree(tfs->in_stream, &node, err)) &&
           node) {
      if (gt_feature_node_try_cast(node) &&
          gt_feature_node_get_attribute(node, "Target")) {
        filter_targetbest((GtFeatureNode*) node, tfs->trees,
                          tfs->target_to_elem);
      }
      else
        gt_dlist_add(tfs->trees, node);
    }
    tfs->next = gt_dlist_first(tfs->trees);
    tfs->in_stream_processed = true;
  }

  if (!had_err) {
    assert(tfs->in_stream_processed);
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

static void targetbest_filter_stream_free(GenomeStream *gs)
{
  TargetbestFilterStream *tfs = targetbest_filter_stream_cast(gs);
  for (; tfs->next != NULL; tfs->next = gt_dlistelem_next(tfs->next))
    gt_genome_node_rec_delete(gt_dlistelem_get_data(tfs->next));
  gt_dlist_delete(tfs->trees);
  hashmap_delete(tfs->target_to_elem);
  genome_stream_delete(tfs->in_stream);
}

const GenomeStreamClass* targetbest_filter_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (TargetbestFilterStream),
                                         targetbest_filter_stream_next_tree,
                                         targetbest_filter_stream_free };
  return &gsc;
}

GenomeStream* targetbest_filter_stream_new(GenomeStream *in_stream)
{
  TargetbestFilterStream *tfs;
  GenomeStream *gs;
  assert(in_stream);
  gs = genome_stream_create(targetbest_filter_stream_class(),
                            genome_stream_is_sorted(in_stream));
  tfs = targetbest_filter_stream_cast(gs);
  tfs->in_stream = genome_stream_ref(in_stream);
  tfs->in_stream_processed = false;
  tfs->trees = gt_dlist_new(NULL);
  tfs->target_to_elem = hashmap_new(HASH_STRING, gt_free_func, NULL);
  return gs;
}
