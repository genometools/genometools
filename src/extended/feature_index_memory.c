/*
  Copyright (c) 2007-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007      Malte Mader <mader@zbh.uni-hamburg.de>
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2007-2012 Center for Bioinformatics, University of Hamburg

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
#include "core/class_alloc_lock.h"
#include "core/cstr_table.h"
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/hashmap.h"
#include "core/interval_tree.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_index_memory.h"
#include "extended/feature_index_rep.h"
#include "extended/feature_index.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"

struct GtFeatureIndexMemory {
  const GtFeatureIndex parent_instance;
  GtHashmap *regions;
  GtHashmap *nodes_in_index;
  GtArray *ids;
  char *firstseqid;
  unsigned long nof_region_nodes,
                reference_count,
                nof_nodes;
};

#define gt_feature_index_memory_cast(FI)\
        gt_feature_index_cast(gt_feature_index_memory_class(), FI)

typedef struct {
  GtIntervalTree *features;
  GtRegionNode *region;
  GtRange dyn_range;
} RegionInfo;

static void region_info_delete(RegionInfo *info)
{
  gt_interval_tree_delete(info->features);
  if (info->region)
    gt_genome_node_delete((GtGenomeNode*)info->region);
  gt_free(info);
}

int gt_feature_index_memory_add_region_node(GtFeatureIndex *gfi,
                                            GtRegionNode *rn,
                                            GT_UNUSED GtError *err)
{
  char *seqid;
  GtFeatureIndexMemory *fi;
  RegionInfo *info;
  fi = gt_feature_index_memory_cast(gfi);
  gt_assert(fi && rn);
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) rn));
  if (!gt_hashmap_get(fi->regions, seqid)) {
    info = gt_calloc(1, sizeof (RegionInfo));
    info->region = (GtRegionNode*) gt_genome_node_ref((GtGenomeNode*) rn);
    info->features = gt_interval_tree_new((GtFree)
                                          gt_genome_node_delete);
    info->dyn_range.start = ~0UL;
    info->dyn_range.end   = 0;
    gt_hashmap_add(fi->regions, seqid, info);
    if (fi->nof_region_nodes++ == 0)
      fi->firstseqid = seqid;
  }
  return 0;
}

int gt_feature_index_memory_add_feature_node(GtFeatureIndex *gfi,
                                             GtFeatureNode *fn,
                                             GT_UNUSED GtError *err)
{
  GtGenomeNode *gn;
  char* seqid;
  GtFeatureIndexMemory *fi;
  GtRange node_range;
  RegionInfo *info;
  GtIntervalTreeNode *new_node;
  gt_assert(gfi && fn);

  fi = gt_feature_index_memory_cast(gfi);
  gn = gt_genome_node_ref((GtGenomeNode*) fn);
  /* assign id number as 'primary key' */
  gt_hashmap_add(fi->nodes_in_index, gn, gn);
  /* get information about seqid and range */
  node_range = gt_genome_node_get_range(gn);
  seqid = gt_str_get(gt_genome_node_get_seqid(gn));
  info = (RegionInfo*) gt_hashmap_get(fi->regions, seqid);

  /* If the seqid was encountered for the first time, no sequence
     region nodes have been visited before. We therefore must create a new
     index entry and maintain our own GtRange. */
  if (!info)
  {
    info = gt_calloc(1, sizeof (RegionInfo));
    info->region = NULL;
    info->features = gt_interval_tree_new((GtFree)
                                          gt_genome_node_delete);
    info->dyn_range.start = ~0UL;
    info->dyn_range.end   = 0;
    gt_hashmap_add(fi->regions, seqid, info);
    if (fi->nof_region_nodes++ == 0)
      fi->firstseqid = seqid;
  }

  /* add node to the appropriate array in the hashtable */
  new_node = gt_interval_tree_node_new(gn, node_range.start, node_range.end);
  gt_interval_tree_insert(info->features, new_node);
  /* update dynamic range */
  info->dyn_range.start = MIN(info->dyn_range.start, node_range.start);
  info->dyn_range.end = MAX(info->dyn_range.end, node_range.end);
  return 0;
}

typedef struct {
  GtIntervalTreeNode *node;
  GtGenomeNode *genome_node;
} GtFeatureIndexMemoryByPtrExtractInfo;

static int gt_feature_index_memory_get_itreenode_by_ptr(GtIntervalTreeNode *n,
                                                        void *data)
{
  GtFeatureIndexMemoryByPtrExtractInfo *i =
                                   (GtFeatureIndexMemoryByPtrExtractInfo*) data;
  if (i->genome_node == gt_interval_tree_node_get_data(n)) {
    i->node = n;
  }
  return 0;
}

int gt_feature_index_memory_remove_node(GtFeatureIndex *gfi,
                                        GtFeatureNode *gn,
                                        GT_UNUSED GtError *err)
{
  char* seqid;
  GtFeatureIndexMemory *fi;
  GtRange node_range;
  GtFeatureIndexMemoryByPtrExtractInfo info;
  RegionInfo *rinfo;
  gt_assert(gfi && gn);

  fi = gt_feature_index_memory_cast(gfi);
  node_range = gt_genome_node_get_range((GtGenomeNode*) gn);
  if (!gt_hashmap_get(fi->nodes_in_index, gn))
    return 0;
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) gn));
  rinfo = (RegionInfo*) gt_hashmap_get(fi->regions, seqid);
  if (!rinfo)
    return 0;
  info.genome_node = (GtGenomeNode*) gn;
  info.node = NULL;

  gt_interval_tree_iterate_overlapping(rinfo->features,
                                   gt_feature_index_memory_get_itreenode_by_ptr,
                                   node_range.start,
                                   node_range.end,
                                   &info);

  if (info.node)
    gt_interval_tree_remove(rinfo->features, info.node);
  return 0;
}

static int collect_features_from_itree(GtIntervalTreeNode *node, void *data)
{
  GtArray *a = (GtArray*) data;
  GtGenomeNode *gn = (GtGenomeNode*) gt_interval_tree_node_get_data(node);
  gt_array_add(a, gn);
  return 0;
}

GtArray* gt_feature_index_memory_get_features_for_seqid(GtFeatureIndex *gfi,
                                                        const char *seqid,
                                                        GT_UNUSED GtError *err)
{
  RegionInfo *ri;
  GT_UNUSED int had_err = 0;
  GtArray *a;
  GtFeatureIndexMemory *fi;
  gt_assert(gfi && seqid);
  fi = gt_feature_index_memory_cast(gfi);
  a = gt_array_new(sizeof (GtFeatureNode*));
  ri = (RegionInfo*) gt_hashmap_get(fi->regions, seqid);
  if (ri) {
    had_err = gt_interval_tree_traverse(ri->features,
                                        collect_features_from_itree,
                                        a);
  }
  gt_assert(!had_err);   /* collect_features_from_itree() is sane */
  return a;
}

static int gt_genome_node_cmp_range_start(const void *v1, const void *v2)
{
  GtGenomeNode *n1, *n2;
  n1 = *(GtGenomeNode**) v1;
  n2 = *(GtGenomeNode**) v2;
  return gt_genome_node_compare(&n1, &n2);
}

int gt_feature_index_memory_get_features_for_range(GtFeatureIndex *gfi,
                                                   GtArray *results,
                                                   const char *seqid,
                                                   const GtRange *qry_range,
                                                   GtError *err)
{
  RegionInfo *ri;
  GtFeatureIndexMemory *fi;
  gt_error_check(err);
  gt_assert(gfi && results);

  fi = gt_feature_index_memory_cast(gfi);
  ri = (RegionInfo*) gt_hashmap_get(fi->regions, seqid);
  if (!ri) {
    gt_error_set(err, "feature index does not contain the given sequence id");
    return -1;
  }
  gt_interval_tree_find_all_overlapping(ri->features, qry_range->start,
                                        qry_range->end, results);
  gt_array_sort(results, gt_genome_node_cmp_range_start);
  return 0;
}

GtFeatureNode*  gt_feature_index_memory_get_node_by_ptr(GtFeatureIndexMemory
                                                                          *fim,
                                                        GtFeatureNode *ptr,
                                                        GtError *err)
{
  GtFeatureNode *retnode;
  gt_assert(fim);

  if (!(retnode = gt_hashmap_get(fim->nodes_in_index, ptr))) {
    gt_error_set(err, "feature index does not contain a node with address %p",
                 ptr);
  }
  return retnode;
}

char* gt_feature_index_memory_get_first_seqid(const GtFeatureIndex *gfi,
                                              GT_UNUSED GtError *err)
{
  GtFeatureIndexMemory *fi;
  gt_assert(gfi);

  fi = gt_feature_index_memory_cast((GtFeatureIndex*) gfi);
  return fi->firstseqid ? gt_cstr_dup(fi->firstseqid) : NULL;
}

static int store_seqid(void *key, GT_UNUSED void *value, void *data,
                       GT_UNUSED GtError *err)
{
  GtCstrTable *seqids = (GtCstrTable*) data;
  const char *seqid = (const char*) key;
  gt_assert(seqids && seqid);
  if (!gt_cstr_table_get(seqids, seqid)) {
    gt_cstr_table_add(seqids, seqid);
  }
  return 0;
}

GtStrArray* gt_feature_index_memory_get_seqids(const GtFeatureIndex *gfi,
                                               GT_UNUSED GtError *err)
{
  GtCstrTable* seqids;
  GtStrArray *ret;
  GT_UNUSED int rval;
  GtFeatureIndexMemory *fi;
  gt_assert(gfi);

  fi = gt_feature_index_memory_cast((GtFeatureIndex*) gfi);
  seqids = gt_cstr_table_new();
  rval = gt_hashmap_foreach_in_key_order(fi->regions, store_seqid, seqids,
                                         NULL);
  gt_assert(!rval); /* store_seqid() is sane */
  ret = gt_cstr_table_get_all(seqids);
  gt_cstr_table_delete(seqids);
  return ret;
}

int gt_feature_index_memory_get_range_for_seqid(GtFeatureIndex *gfi,
                                                GtRange *range,
                                                const char *seqid,
                                                GT_UNUSED GtError *err)
{
  RegionInfo *info;
  GtFeatureIndexMemory *fi;
  gt_assert(gfi && range && seqid);
  fi = gt_feature_index_memory_cast(gfi);
  info = (RegionInfo*) gt_hashmap_get(fi->regions, seqid);
  gt_assert(info);

  if (info->dyn_range.start != ~0UL && info->dyn_range.end != 0) {
    range->start = info->dyn_range.start;
    range->end = info->dyn_range.end;
  }
  else if (info->region)
    *range = gt_genome_node_get_range((GtGenomeNode*) info->region);
  return 0;
}

int gt_feature_index_memory_has_seqid(const GtFeatureIndex *gfi,
                                      bool *has_seqid,
                                      const char *seqid,
                                      GT_UNUSED GtError *err)
{
  GtFeatureIndexMemory *fi;
  gt_assert(gfi);

  fi = gt_feature_index_memory_cast((GtFeatureIndex*) gfi);
  *has_seqid = (gt_hashmap_get(fi->regions, seqid) != NULL);
  return 0;
}

void gt_feature_index_memory_delete(GtFeatureIndex *gfi)
{
  GtFeatureIndexMemory *fi;
  if (!gfi) return;
  fi = gt_feature_index_memory_cast(gfi);
  gt_hashmap_delete(fi->regions);
  gt_hashmap_delete(fi->nodes_in_index);
}

const GtFeatureIndexClass* gt_feature_index_memory_class(void)
{
  static const GtFeatureIndexClass *fic = NULL;
  gt_class_alloc_lock_enter();
  if (!fic) {
    fic = gt_feature_index_class_new(sizeof (GtFeatureIndexMemory),
                     gt_feature_index_memory_add_region_node,
                     gt_feature_index_memory_add_feature_node,
                     gt_feature_index_memory_remove_node,
                     gt_feature_index_memory_get_features_for_seqid,
                     gt_feature_index_memory_get_features_for_range,
                     gt_feature_index_memory_get_first_seqid,
                     NULL,
                     gt_feature_index_memory_get_seqids,
                     gt_feature_index_memory_get_range_for_seqid,
                     gt_feature_index_memory_has_seqid,
                     gt_feature_index_memory_delete);
  }
  gt_class_alloc_lock_leave();
  return fic;
}

GtFeatureIndex* gt_feature_index_memory_new(void)
{
  GtFeatureIndexMemory *fim;
  GtFeatureIndex *fi;
  fi = gt_feature_index_create(gt_feature_index_memory_class());
  fim = gt_feature_index_memory_cast(fi);
  fim->nof_nodes = 0;
  fim->regions = gt_hashmap_new(GT_HASH_STRING, NULL,
                                (GtFree) region_info_delete);
  fim->nodes_in_index = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  return fi;
}

int gt_feature_index_memory_unit_test(GtError *err)
{
  int had_err = 0, status = 0;
  GtFeatureIndex *fi = NULL;
  GtFeatureNode *tmp, *fn;
  GtError *testerr;
  gt_error_check(err);

  fi = gt_feature_index_memory_new();
  gt_ensure(had_err, fi);

  /* run generic feature index tests */
  had_err = gt_feature_index_unit_test(fi, err);
  gt_ensure(had_err, status == 0);

  /* run subclass specific tests */
  testerr = gt_error_new();
  fn = gt_feature_node_cast(gt_feature_node_new_standard_gene());
  gt_ensure(had_err, !gt_feature_index_add_feature_node(fi, fn, testerr));

  /* test gt_feature_index_memory_get_node_by_ptr() */
  tmp = gt_feature_index_memory_get_node_by_ptr(
                                            gt_feature_index_memory_cast(fi),
                                            fn, testerr);
  gt_ensure(had_err, tmp == fn);
  gt_ensure(had_err, !gt_error_is_set(testerr));

  /* we don't store NULL pointers */
  tmp = gt_feature_index_memory_get_node_by_ptr(
                                               gt_feature_index_memory_cast(fi),
                                               (GtFeatureNode*) 0,
                                                testerr);
  gt_ensure(had_err, tmp == NULL);
  gt_ensure(had_err, gt_error_is_set(testerr));
  gt_genome_node_delete((GtGenomeNode*) fn);
  gt_feature_index_delete(fi);

  gt_error_delete(testerr);
  return had_err;
}
