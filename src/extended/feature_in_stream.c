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
  GtStrArray *seqids;
  GtQueue *regioncache;
  GtArray *featurecache;
  GtUword seqindex;
  bool useorig;
  bool init;
};

#define feature_in_stream_cast(GS)\
        gt_node_stream_cast(gt_feature_in_stream_class(), GS)

void feature_in_stream_init(GtFeatureInStream *stream)
{
  GtUword i;
  GtError *error = gt_error_new();

  stream->seqids = gt_feature_index_get_seqids(stream->fi, error);
  stream->seqindex = 0;
  for (i = 0; i < gt_str_array_size(stream->seqids); i++)
  {
    const char *seqid = gt_str_array_get(stream->seqids, i);
    GtRange seqrange;
    if (stream->useorig)
      gt_feature_index_get_orig_range_for_seqid(stream->fi, &seqrange, seqid,
                                                error);
    else
      gt_feature_index_get_range_for_seqid(stream->fi, &seqrange, seqid, error);
    GtStr *seqstr = gt_str_new_cstr(seqid);
    GtGenomeNode *rn = gt_region_node_new(seqstr, seqrange.start, seqrange.end);
    gt_queue_add(stream->regioncache, rn);
    gt_str_delete(seqstr);
  }
  gt_error_delete(error);
}

static int feature_in_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                   GtError *error)
{
  GtFeatureInStream *stream = feature_in_stream_cast(ns);
  gt_error_check(error);

  if (!stream->init)
  {
    feature_in_stream_init(stream);
    stream->init = true;
  }

  if (gt_queue_size(stream->regioncache) > 0)
  {
    GtGenomeNode *region = gt_queue_get(stream->regioncache);
    *gn = region;
    return 0;
  }

  if (stream->featurecache == NULL || gt_array_size(stream->featurecache) == 0)
  {
    GtUword numfeats = 0;
    if (stream->featurecache != NULL)
    {
      gt_array_delete(stream->featurecache);
      stream->featurecache = NULL;
    }

    if (stream->seqindex == gt_str_array_size(stream->seqids))
    {
      *gn = NULL;
      return 0;
    }

    do
    {
      const char *seqid = gt_str_array_get(stream->seqids, stream->seqindex++);
      stream->featurecache = gt_feature_index_get_features_for_seqid(stream->fi,
                                                                     seqid,
                                                                     error);
      numfeats = gt_array_size(stream->featurecache);
    } while (numfeats == 0 &&
             stream->seqindex < gt_str_array_size(stream->seqids));
    if (numfeats > 0)
    {
      gt_array_sort(stream->featurecache, (GtCompare)gt_genome_node_compare);
      gt_array_reverse(stream->featurecache);
    }
  }

  if (gt_array_size(stream->featurecache) == 0)
  {
    *gn = NULL;
  }
  else
  {
    GtGenomeNode *feat = *(GtGenomeNode **)gt_array_pop(stream->featurecache);
    *gn = gt_genome_node_ref(feat);
  }
  return 0;
}

static void feature_in_stream_free(GtNodeStream *ns)
{
  GtFeatureInStream *stream = feature_in_stream_cast(ns);
  gt_str_array_delete(stream->seqids);
  while (gt_queue_size(stream->regioncache) > 0)
  {
    GtGenomeNode *gn = gt_queue_get(stream->regioncache);
    gt_genome_node_delete(gn);
  }
  gt_queue_delete(stream->regioncache);
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
  stream->regioncache = gt_queue_new();
  stream->featurecache = NULL;
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

  return fi;
}

int gt_feature_in_stream_unit_test(GtError *error)
{
  GtNodeStream *src, *dest;
  GtFeatureIndex *prefeat, *postfeat;
  GtRange range1, range1test, range2, range2test;

  prefeat = in_stream_test_data(error);
  postfeat = gt_feature_index_memory_new();
  src = gt_feature_in_stream_new(prefeat);
  dest = gt_feature_out_stream_new(src, postfeat);
  int result = gt_node_stream_pull(dest, error);
  if (result == -1)
    return -1;

  GtStrArray *seqids = gt_feature_index_get_seqids(postfeat, error);
  if (gt_str_array_size(seqids) != 2)
  {
    gt_error_set(error, "error in feature_in_stream unit test 1: expected 2 "
                 "seqids, found "GT_WU"", gt_str_array_size(seqids));
    return -1;
  }
  gt_str_array_delete(seqids);

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

  return 0;
}
