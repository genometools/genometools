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
#include "extended/feature_index.h"
#include "extended/feature_in_stream.h"
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
};

#define feature_in_stream_cast(GS)\
        gt_node_stream_cast(gt_feature_in_stream_class(), GS)

static int feature_in_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                   GtError *error)
{
  GtFeatureInStream *stream = feature_in_stream_cast(ns);
  gt_error_check(error);

  if (gt_queue_size(stream->regioncache) > 0)
  {
    GtGenomeNode *region = gt_queue_get(stream->regioncache);
    *gn = region;
    return 0;
  }

  if (stream->featurecache == NULL || gt_array_size(stream->featurecache) == 0)
  {
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

    const char *seqid = gt_str_array_get(stream->seqids, stream->seqindex++);
    stream->featurecache = gt_feature_index_get_features_for_seqid(stream->fi,
                                                                   seqid,
                                                                   error);
    gt_array_sort(stream->featurecache, (GtCompare)gt_genome_node_compare);
    gt_array_reverse(stream->featurecache);
  }

  GtGenomeNode *feat = *(GtGenomeNode **)gt_array_pop(stream->featurecache);
  *gn = feat;
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
    gt_feature_index_get_range_for_seqid(stream->fi, &seqrange, seqid, error);
    GtStr *seqstr = gt_str_new_cstr(seqid);
    GtGenomeNode *rn = gt_region_node_new(seqstr, seqrange.start, seqrange.end);
    gt_queue_add(stream->regioncache, rn);
    gt_str_delete(seqstr);
  }
  gt_error_delete(error);
}

GtNodeStream* gt_feature_in_stream_new(GtFeatureIndex *fi)
{
  GtNodeStream *ns;
  GtFeatureInStream *stream;
  ns = gt_node_stream_create(gt_feature_in_stream_class(), true);
  stream = feature_in_stream_cast(ns);
  stream->fi = fi;
  stream->regioncache = gt_queue_new();
  stream->featurecache = NULL;
  feature_in_stream_init(stream);
  return ns;
}
