/*
  Copyright (c) 2007-2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007      Malte Mader <mmader@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007      Chr. Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
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

#include <string.h>
#include "annotationsketch/feature_index.h"
#include "annotationsketch/feature_stream.h"
#include "annotationsketch/feature_visitor.h"
#include "core/ensure.h"
#include "core/hashmap.h"
#include "core/interval_tree.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/queue.h"
#include "core/range.h"
#include "core/undef.h"
#include "core/unused.h"
#include "extended/feature_type_factory_builtin.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"

struct GT_FeatureIndex {
  Hashmap *regions;
  char *firstseqid;
  unsigned int nof_sequence_regions,
               reference_count;
};

typedef struct {
  IntervalTree *features;
  GT_SequenceRegion *region;
  GT_Range dyn_range;
} RegionInfo;

static void region_info_delete(RegionInfo *info)
{
  interval_tree_delete(info->features);
  if (info->region)
    gt_genome_node_delete((GT_GenomeNode*)info->region);
  ma_free(info);
}

GT_FeatureIndex* gt_feature_index_new(void)
{
  GT_FeatureIndex *fi;
  fi = ma_calloc(1, sizeof (GT_FeatureIndex));
  fi->regions = hashmap_new(HASH_STRING, NULL, (FreeFunc) region_info_delete);
  return fi;
}

GT_FeatureIndex* gt_feature_index_ref(GT_FeatureIndex *fi)
{
  assert(fi);
  fi->reference_count++;
  return fi;
}

void gt_feature_index_add_sequence_region(GT_FeatureIndex *fi,
                                          GT_SequenceRegion *sr)
{
  char *seqid;
  RegionInfo *info;
  assert(fi && sr);
  seqid = str_get(gt_genome_node_get_seqid((GT_GenomeNode*) sr));
  if (!hashmap_get(fi->regions, seqid)) {
    info = ma_malloc(sizeof (RegionInfo));
    info->region = (GT_SequenceRegion*) gt_genome_node_ref((GT_GenomeNode*) sr);
    info->features = interval_tree_new((FreeFunc) gt_genome_node_rec_delete);
    info->dyn_range.start = ~0UL;
    info->dyn_range.end   = 0;
    hashmap_add(fi->regions, seqid, info);
    if (fi->nof_sequence_regions++ == 0)
      fi->firstseqid = seqid;
  }
}

void gt_feature_index_add_genome_feature(GT_FeatureIndex *fi, GT_GenomeFeature *gf)
{
  GT_GenomeNode *gn;
  char* seqid;
  GT_Range node_range;
  RegionInfo *info;

  assert(fi && gf);

  gn = gt_genome_node_rec_ref((GT_GenomeNode*) gf);
  /* get information about seqid and range */
  node_range = gt_genome_node_get_range(gn);
  seqid = str_get(gt_genome_node_get_seqid(gn));
  info = (RegionInfo*) hashmap_get(fi->regions, seqid);

  /* If the seqid was encountered for the first time, no sequence
     region nodes have been visited before. We therefore must create a new
     index entry and maintain our own GT_Range. */
  if (!info)
  {
    info = ma_calloc(1, sizeof (RegionInfo));
    info->region = NULL;
    info->features = interval_tree_new((FreeFunc) gt_genome_node_rec_delete);
    info->dyn_range.start = ~0UL;
    info->dyn_range.end   = 0;
    hashmap_add(fi->regions, seqid, info);
    if (fi->nof_sequence_regions++ == 0)
      fi->firstseqid = seqid;
  }

  /* add node to the appropriate array in the hashtable */
  IntervalTreeNode *new_node = interval_tree_node_new(gn, node_range.start,
                                                      node_range.end);
  interval_tree_insert(info->features, new_node);
  /* update dynamic range */
  info->dyn_range.start = MIN(info->dyn_range.start, node_range.start);
  info->dyn_range.end = MAX(info->dyn_range.end, node_range.end);
}

int gt_feature_index_add_gff3file(GT_FeatureIndex *feature_index,
                               const char *gff3file, GT_Error *err)
{
  GenomeStream *gff3_in_stream;
  GT_GenomeNode *gn;
  Queue *queue;
  int had_err = 0;
  gt_error_check(err);
  assert(feature_index && gff3file);
  queue = queue_new();
  gff3_in_stream = gff3_in_stream_new_unsorted(1, &gff3file, false, false);
  while (!(had_err = genome_stream_next_tree(gff3_in_stream, &gn, err)) && gn)
    queue_add(queue, gn);
  if (!had_err) {
    GenomeVisitor *feature_visitor = feature_visitor_new(feature_index);
    while (queue_size(queue)) {
      gn = queue_get(queue);
      had_err = gt_genome_node_accept(gn, feature_visitor, NULL);
      assert(!had_err); /* cannot happen */
    }
    genome_visitor_delete(feature_visitor);
  }
  genome_stream_delete(gff3_in_stream);
  while (queue_size(queue))
    gt_genome_node_rec_delete(queue_get(queue));
  queue_delete(queue);
  return had_err;
}

static int collect_features_from_itree(IntervalTreeNode *node, void *data)
{
  GT_Array *a = (GT_Array*) data;
  GT_GenomeNode *gn = (GT_GenomeNode*) interval_tree_node_get_data(node);
  gt_array_add(a, gn);
  return 0;
}

GT_Array* gt_feature_index_get_features_for_seqid(GT_FeatureIndex *fi,
                                                  const char *seqid)
{
  RegionInfo *ri;
  int had_err = 0;
  GT_Array *a;
  assert(fi && seqid);
  a = gt_array_new(sizeof (GT_GenomeFeature*));
  ri = (RegionInfo*) hashmap_get(fi->regions, seqid);
  if (ri)
    had_err = interval_tree_traverse(ri->features,
                                     collect_features_from_itree,
                                     a);
  assert(!had_err);   /* collect_features_from_itree() is sane */
  return a;
}

static int gt_genome_node_cmp_range_start(const void *v1, const void *v2)
{
  GT_GenomeNode *n1, *n2;
  n1 = *(GT_GenomeNode**) v1;
  n2 = *(GT_GenomeNode**) v2;
  return gt_genome_node_compare(&n1, &n2);
}

int gt_feature_index_get_features_for_range(GT_FeatureIndex *fi, GT_Array *results,
                                         const char *seqid, GT_Range qry_range,
                                         GT_Error *err)
{
  RegionInfo *ri;
  gt_error_check(err);
  assert(fi && results);

  ri = (RegionInfo*) hashmap_get(fi->regions, seqid);
  if (!ri) {
    gt_error_set(err, "feature index does not contain the given sequence id");
    return -1;
  }
  interval_tree_find_all_overlapping(ri->features, qry_range.start,
                                     qry_range.end, results);
  gt_array_sort(results, gt_genome_node_cmp_range_start);
  return 0;
}

const char* gt_feature_index_get_first_seqid(const GT_FeatureIndex *fi)
{
  assert(fi);
  return fi->firstseqid;
}

int store_seqid(void *key, UNUSED void *value, void *data, UNUSED GT_Error *err)
{
  GT_StrArray *seqids = (GT_StrArray*) data;
  const char *seqid = (const char*) key;
  assert(seqids && seqid);
  gt_strarray_add_cstr(seqids, seqid);
  return 0;
}

GT_StrArray* gt_feature_index_get_seqids(const GT_FeatureIndex *fi)
{
  GT_StrArray* seqids;
  int rval;
  assert(fi);
  seqids = gt_strarray_new();
  rval = hashmap_foreach_in_key_order(fi->regions, store_seqid, seqids, NULL);
  assert(!rval); /* store_seqid() is sane */
  return seqids;
}

void gt_feature_index_get_range_for_seqid(GT_FeatureIndex *fi, GT_Range *range,
                                          const char *seqid)
{
  RegionInfo *info;
  assert(fi && range && seqid);
  info = (RegionInfo*) hashmap_get(fi->regions, seqid);
  assert(info);

  if (info->dyn_range.start != ~0UL && info->dyn_range.end != 0) {
    range->start = info->dyn_range.start;
    range->end = info->dyn_range.end;
  }
  else if (info->region)
    *range = gt_genome_node_get_range((GT_GenomeNode*) info->region);
}

bool gt_feature_index_has_seqid(const GT_FeatureIndex *fi, const char *seqid)
{
  assert(fi);
  return (hashmap_get(fi->regions, seqid));
}

int gt_feature_index_unit_test(GT_Error *err)
{
  FeatureTypeFactory *feature_type_factory;
  GT_GenomeFeatureType *type;
  GT_GenomeNode *gn1, *gn2, *ex1, *ex2, *ex3, *cds1;
  GT_FeatureIndex *fi;
  GT_Range r1, r2, r3, r4, r5, check_range, rs;
  Str *seqid1, *seqid2;
  GT_StrArray *seqids = NULL;
  GT_SequenceRegion *sr1, *sr2;
  GT_Array *features = NULL;
  int had_err = 0;
  gt_error_check(err);

  feature_type_factory = feature_type_factory_builtin_new();

  /* generating some ranges */
  r1.start=100UL; r1.end=1000UL;
  r2.start=100UL; r2.end=300UL;
  r3.start=500UL; r3.end=1000UL;
  r4.start=600UL; r4.end=1200UL;
  r5.start=600UL; r5.end=1000UL;
  rs.start=100UL; rs.end=1200UL;

  /* generating sequnce ids as C-strings */
  seqid1 = str_new_cstr("test1");
  seqid2 = str_new_cstr("test2");

  sr1 = (GT_SequenceRegion*) gt_sequence_regionnew(seqid1, rs);
  sr2 = (GT_SequenceRegion*) gt_sequence_regionnew(seqid2, rs);

  /* generate a new genome feature */
  type = feature_type_factory_create_gft(feature_type_factory, gft_gene);
  gn1 = gt_genome_feature_new(seqid1, type, r1, STRAND_UNKNOWN);

  gn2 = gt_genome_feature_new(seqid2, type, r4, STRAND_UNKNOWN);

  type = feature_type_factory_create_gft(feature_type_factory, gft_exon);
  ex1 = gt_genome_feature_new(seqid1, type, r2, STRAND_UNKNOWN);

  ex2 = gt_genome_feature_new(seqid1, type, r3, STRAND_UNKNOWN);

  ex3 = gt_genome_feature_new(seqid2, type, r4, STRAND_UNKNOWN);

  type = feature_type_factory_create_gft(feature_type_factory, gft_CDS);
  cds1 = gt_genome_feature_new(seqid2, type, r5, STRAND_UNKNOWN);

  /* Determine the structure of our feature tree */
  gt_genome_node_is_part_of_genome_node(gn1, ex1);
  gt_genome_node_is_part_of_genome_node(gn1, ex2);
  gt_genome_node_is_part_of_genome_node(gn2, ex3);
  gt_genome_node_is_part_of_genome_node(gn2, cds1);

  /* create a new feature index on which we can perform some tests */
  fi = gt_feature_index_new();

  ensure(had_err, fi);
  ensure(had_err, !gt_feature_index_has_seqid(fi, "test1"));
  ensure(had_err, !gt_feature_index_has_seqid(fi, "test2"));

  /* add a sequence region directly and check if it has been added */
  gt_feature_index_add_sequence_region(fi, sr1);
  ensure(had_err, gt_feature_index_has_seqid(fi, "test1"));
  ensure(had_err, !gt_feature_index_has_seqid(fi, "test2"));

  gt_feature_index_get_range_for_seqid(fi, &check_range, "test1");
  ensure(had_err, check_range.start == 100UL && check_range.end == 1200UL);

  /* tests if we get a empty data structure for every added sequence region*/
  if (!had_err)
    features = gt_feature_index_get_features_for_seqid(fi, "test1");
  ensure(had_err, features);
  ensure(had_err, gt_array_size(features) == 0);
  gt_array_delete(features);
  features = NULL;

  if (!had_err)
    features = gt_feature_index_get_features_for_seqid(fi, "test2");
  ensure(had_err, features);
  ensure(had_err, gt_array_size(features) == 0);
  gt_array_delete(features);
  features = NULL;

  /* add features to every sequence region and test if the according
     datastructures are not empty anymore. As we have added one genome_feature
     to every sequence region the size has to be one. */
  if (!had_err) {
    gt_feature_index_add_genome_feature(fi, (GT_GenomeFeature*) gn1);
    features = gt_feature_index_get_features_for_seqid(fi, "test1");
  }
  ensure(had_err, gt_array_size(features) == 1UL);
  gt_array_delete(features);
  features = NULL;

  if (!had_err) {
    gt_feature_index_add_genome_feature(fi, (GT_GenomeFeature*) gn2);
    features = gt_feature_index_get_features_for_seqid(fi, "test2");
  }
  ensure(had_err, gt_array_size(features) == 1UL);
  gt_array_delete(features);
  features = NULL;

  /* test gt_feature_index_get_first_seqid() */
  ensure(had_err, gt_feature_index_get_first_seqid(fi));
  ensure(had_err, strcmp("test1", gt_feature_index_get_first_seqid(fi)) == 0);

  if (!had_err) {
    seqids = gt_feature_index_get_seqids(fi);
    ensure(had_err, gt_strarray_size(seqids) == 2);
    ensure(had_err, !strcmp(gt_strarray_get(seqids, 0), "test1"));
    ensure(had_err, !strcmp(gt_strarray_get(seqids, 1), "test2"));
  }

  gt_feature_index_get_range_for_seqid(fi, &check_range, "test1");
  ensure(had_err, check_range.start == 100UL && check_range.end == 1000UL);

  gt_feature_index_get_range_for_seqid(fi, &check_range, "test2");
  ensure(had_err, check_range.start == 600UL && check_range.end == 1200UL);

  if (!had_err)
    features = gt_feature_index_get_features_for_seqid(fi, "test1");
  ensure(had_err, features);
  gt_array_delete(features);

  /* delete all generated objects */
  gt_strarray_delete(seqids);
  gt_feature_index_delete(fi);
  gt_genome_node_rec_delete(gn1);
  gt_genome_node_rec_delete(gn2);
  gt_genome_node_rec_delete((GT_GenomeNode*) sr1);
  gt_genome_node_rec_delete((GT_GenomeNode*) sr2);
  str_delete(seqid1);
  str_delete(seqid2);
  feature_type_factory_delete(feature_type_factory);
  return had_err;
}

void gt_feature_index_delete(GT_FeatureIndex *fi)
{
  if (!fi) return;
  if (fi->reference_count) {
    fi->reference_count--;
    return;
  }
  hashmap_delete(fi->regions);
  ma_free(fi);
}
