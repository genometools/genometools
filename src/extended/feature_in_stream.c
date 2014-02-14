/*
  Copyright (c) 2013-2014 Daniel S. Standage <daniel.standage@gmail.com>

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

#include "core/class_alloc_lock.h"
#include "core/queue_api.h"
#include "core/unused_api.h"
#include "extended/feature_index_memory_api.h"
#include "extended/feature_in_stream.h"
#include "extended/feature_out_stream_api.h"
#include "extended/feature_visitor.h"
#include "extended/genome_node.h"
#include "extended/visitor_stream_api.h"

struct GtFeatureInStream
{
  const GtNodeStream parent_instance;
  GtFeatureIndex *fi;
  GtQueue *cache;
  bool useorig;
  bool init;
};

#define feature_in_stream_cast(GS)\
        gt_node_stream_cast(gt_feature_in_stream_class(), GS)

void feature_in_stream_init(GtFeatureInStream *stream)
{
  GtError *error;
  GtStrArray *seqids;
  GtUword i;

  error = gt_error_new();
  seqids = gt_feature_index_get_seqids(stream->fi, error);

  /* Load all region nodes into the cache */
  for (i = 0; i < gt_str_array_size(seqids); i++)
  {
    GtRange seqrange;
    GtGenomeNode *rn;
    GtStr *seqstr;
    const char *seqid;

    seqid = gt_str_array_get(seqids, i);
    if (stream->useorig)
    {
      gt_feature_index_get_orig_range_for_seqid(stream->fi, &seqrange, seqid,
                                                error);
    }
    else
    {
      gt_feature_index_get_range_for_seqid(stream->fi, &seqrange, seqid, error);
    }

    seqstr = gt_str_new_cstr(seqid);
    rn = gt_region_node_new(seqstr, seqrange.start, seqrange.end);
    gt_queue_add(stream->cache, rn);
    gt_str_delete(seqstr);
  }

  /* Load all feature nodes into the cache */
  for (i = 0; i < gt_str_array_size(seqids); i++)
  {
    GtArray *features;
    const char *seqid;

    seqid = gt_str_array_get(seqids, i);
    features = gt_feature_index_get_features_for_seqid(stream->fi, seqid,error);
    if (gt_array_size(features) > 0)
    {
      gt_array_sort(features, (GtCompare)gt_genome_node_compare);
      gt_array_reverse(features);
      while (gt_array_size(features) > 0)
      {
        GtGenomeNode **fn = gt_array_pop(features);
        gt_queue_add(stream->cache, *fn);
      }
    }
    gt_array_delete(features);
  }

  gt_error_delete(error);
  gt_str_array_delete(seqids);
}

static int feature_in_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                  GT_UNUSED GtError *error)
{
  GtFeatureInStream *stream = feature_in_stream_cast(ns);
  gt_error_check(error);

  if (!stream->init)
  {
    feature_in_stream_init(stream);
    stream->init = true;
  }

  if (gt_queue_size(stream->cache) > 0)
  {
    GtGenomeNode *node = gt_queue_get(stream->cache);
    GtFeatureNode *fn = gt_feature_node_try_cast(node);
    if (fn)
      gt_genome_node_ref(node);
    *gn = node;
    return 0;
  }

  *gn = NULL;
  return 0;
}

static void feature_in_stream_free(GtNodeStream *ns)
{
  GtFeatureInStream *stream = feature_in_stream_cast(ns);
  while (gt_queue_size(stream->cache) > 0)
  {
    GtGenomeNode *gn = gt_queue_get(stream->cache);
    gt_genome_node_delete(gn);
  }
  gt_queue_delete(stream->cache);
}

const GtNodeStreamClass *gt_feature_in_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (GtFeatureInStream),
                                   feature_in_stream_free,
                                   feature_in_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

void gt_feature_in_stream_use_orig_ranges(GtFeatureInStream *stream)
{
  stream->useorig = true;
}

GtNodeStream* gt_feature_in_stream_new(GtFeatureIndex *fi)
{
  GtNodeStream *ns;
  GtFeatureInStream *stream;
  gt_assert(fi);

  ns = gt_node_stream_create(gt_feature_in_stream_class(), true);
  stream = feature_in_stream_cast(ns);
  stream->fi = fi;
  stream->cache = gt_queue_new();
  stream->useorig = false;
  stream->init = false;
  return ns;
}

static GtFeatureIndex *in_stream_test_data(GtError *error)
{
  GtFeatureIndex *fi;
  GtFeatureNode *fn;
  GtGenomeNode *gn;
  GtRegionNode *rn;
  GtStr *seqid;

  fi = gt_feature_index_memory_new();

  seqid = gt_str_new_cstr("chr1");
  gn = gt_region_node_new(seqid, 1, 100000);
  rn = gt_region_node_cast(gn);
  gt_feature_index_add_region_node(fi, rn, error);
  gt_genome_node_delete(gn);

  gn = gt_feature_node_new(seqid, "region", 500, 5000, GT_STRAND_BOTH);
  fn = gt_feature_node_cast(gn);
  gt_feature_index_add_feature_node(fi, fn, error);
  gt_genome_node_delete(gn);

  gn = gt_feature_node_new(seqid, "region", 50000, 75000, GT_STRAND_BOTH);
  fn = gt_feature_node_cast(gn);
  gt_feature_index_add_feature_node(fi, fn, error);
  gt_genome_node_delete(gn);

  gt_str_delete(seqid);
  seqid = gt_str_new_cstr("scf0001");
  gn = gt_region_node_new(seqid, 1, 10000);
  rn = gt_region_node_cast(gn);
  gt_feature_index_add_region_node(fi, rn, error);
  gt_genome_node_delete(gn);

  gn = gt_feature_node_new(seqid, "mRNA", 4000, 6000, GT_STRAND_REVERSE);
  fn = gt_feature_node_cast(gn);
  gt_feature_index_add_feature_node(fi, fn, error);
  gt_genome_node_delete(gn);

  gn = gt_feature_node_new(seqid, "mRNA", 7000, 9500, GT_STRAND_FORWARD);
  fn = gt_feature_node_cast(gn);
  gt_feature_index_add_feature_node(fi, fn, error);
  gt_genome_node_delete(gn);

  gt_str_delete(seqid);
  seqid = gt_str_new_cstr("ChrM");
  gn = gt_region_node_new(seqid, 1, 25000);
  rn = gt_region_node_cast(gn);
  gt_feature_index_add_region_node(fi, rn, error);
  gt_genome_node_delete(gn);

  gt_str_delete(seqid);

  return fi;
}

int gt_feature_in_stream_unit_test(GtError *error)
{
  GtNodeStream *src, *dest;
  GtFeatureIndex *prefeat, *postfeat;
  GtRange range1, range1test, range2, range2test, range3, range3test;
  GtStrArray *seqids;

  prefeat = in_stream_test_data(error);
  postfeat = gt_feature_index_memory_new();
  src = gt_feature_in_stream_new(prefeat);
  dest = gt_feature_out_stream_new(src, postfeat);
  int result = gt_node_stream_pull(dest, error);
  if (result == -1)
    return -1;

  seqids = gt_feature_index_get_seqids(postfeat, error);
  if (gt_str_array_size(seqids) != 3)
  {
    gt_error_set(error, "error in feature_in_stream unit test 1: expected 3 "
                 "seqids, found "GT_WU"", gt_str_array_size(seqids));
    return -1;
  }
  gt_str_array_delete(seqids);

  range1test.start = 500;  range1test.end = 75000;
  range2test.start = 4000; range2test.end = 9500;
  gt_feature_index_get_range_for_seqid(postfeat, &range1, "chr1", error);
  gt_feature_index_get_range_for_seqid(postfeat, &range2, "scf0001", error);
  if (gt_range_compare(&range1, &range1test) ||
      gt_range_compare(&range2, &range2test))
  {
    gt_error_set(error, "error in feature_in_stream unit test 1: incorrect "
                 "sequence regions");
    return -1;
  }

  gt_feature_index_get_orig_range_for_seqid(postfeat, &range1, "chr1", error);
  gt_feature_index_get_orig_range_for_seqid(postfeat, &range2, "scf0001",error);
  if (gt_range_compare(&range1, &range1test) ||
      gt_range_compare(&range2, &range2test))
  {
    gt_error_set(error, "error in feature_in_stream unit test 1: incorrect "
                 "sequence regions");
    return -1;
  }
  gt_feature_index_delete(prefeat);
  gt_feature_index_delete(postfeat);
  gt_node_stream_delete(src);
  gt_node_stream_delete(dest);

  prefeat = in_stream_test_data(error);
  postfeat = gt_feature_index_memory_new();
  src = gt_feature_in_stream_new(prefeat);
  dest = gt_feature_out_stream_new(src, postfeat);
  gt_feature_in_stream_use_orig_ranges((GtFeatureInStream *)src);
  result = gt_node_stream_pull(dest, error);
  if (result == -1)
    return -1;

  range1test.start = 500;  range1test.end = 75000;
  range2test.start = 4000; range2test.end = 9500;
  gt_feature_index_get_range_for_seqid(postfeat, &range1, "chr1", error);
  gt_feature_index_get_range_for_seqid(postfeat, &range2, "scf0001",error);
  if (gt_range_compare(&range1, &range1test) ||
      gt_range_compare(&range2, &range2test))
  {
    gt_error_set(error, "error in feature_in_stream unit test 1: incorrect "
                 "sequence regions");
    return -1;
  }
  range1test.start = 1; range1test.end = 100000;
  range2test.start = 1; range2test.end = 10000;
  range3test.start = 1; range3test.end = 25000;
  gt_feature_index_get_orig_range_for_seqid(postfeat, &range1, "chr1", error);
  gt_feature_index_get_orig_range_for_seqid(postfeat, &range2, "scf0001",error);
  gt_feature_index_get_orig_range_for_seqid(postfeat, &range3, "ChrM", error);
  if (gt_range_compare(&range1, &range1test) ||
      gt_range_compare(&range2, &range2test) ||
      gt_range_compare(&range3, &range3test))
  {
    gt_error_set(error, "error in feature_in_stream unit test 1: incorrect "
                 "sequence regions");
    return -1;
  }
  gt_feature_index_delete(prefeat);
  gt_feature_index_delete(postfeat);
  gt_node_stream_delete(src);
  gt_node_stream_delete(dest);

  return 0;
}
