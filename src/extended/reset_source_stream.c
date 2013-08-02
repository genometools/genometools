/*
  Copyright (c) 2013 Daniel Standage <daniel.standage@gmail.com>

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
#include "core/undef_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_stream_api.h"
#include "extended/reset_source_stream.h"

struct GtResetSourceStream{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtStr *newsource;
};

#define reset_source_stream_cast(GS)\
        gt_node_stream_cast(gt_reset_source_stream_class(), GS)

static int reset_source_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                    GtError *err)
{
  GtResetSourceStream *rss;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(err);
  rss = reset_source_stream_cast(ns);

  had_err = gt_node_stream_next(rss->in_stream, gn, err);
  if(had_err)
    return had_err;
  if(!*gn)
    return 0;
  
  fn = gt_feature_node_try_cast(*gn);
  if(!fn)
    return 0;

  GtFeatureNode *current;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    gt_feature_node_set_source(current, rss->newsource);
  }
  gt_feature_node_iterator_delete(iter);

  return 0;
}

static void reset_source_stream_free(GtNodeStream *ns)
{
  GtResetSourceStream *rss = reset_source_stream_cast(ns);
  gt_node_stream_delete(rss->in_stream);
  gt_str_delete(rss->newsource);
}

const GtNodeStreamClass *gt_reset_source_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtResetSourceStream),
                                   reset_source_stream_free,
                                   reset_source_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_reset_source_stream_new(GtNodeStream *in_stream,
                                         GtStr *newsource)
{
  GtNodeStream *ns;
  GtResetSourceStream *rss;
  gt_assert(in_stream && newsource);
  ns = gt_node_stream_create(gt_reset_source_stream_class(), true);
  rss = reset_source_stream_cast(ns);
  rss->in_stream = gt_node_stream_ref(in_stream);
  rss->newsource = gt_str_ref(newsource);
  return ns;
}
