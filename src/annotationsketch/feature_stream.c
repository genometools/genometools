/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "annotationsketch/feature_stream.h"
#include "annotationsketch/feature_visitor.h"
#include "annotationsketch/feature_index.h"
#include "extended/node_stream_api.h"

struct GtFeatureStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *feature_visitor;
};

#define feature_stream_cast(GS)\
        gt_node_stream_cast(gt_feature_stream_class(), GS)

static int feature_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                               GtError *err)
{
  GtFeatureStream *feature_stream;
  int had_err;
  gt_error_check(err);
  feature_stream = feature_stream_cast(gs);
  had_err = gt_node_stream_next(feature_stream->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, feature_stream->feature_visitor, err);
  return had_err;
}

static void feature_stream_free(GtNodeStream *gs)
{
  GtFeatureStream *feature_stream = feature_stream_cast(gs);
  gt_node_stream_delete(feature_stream->in_stream);
  gt_node_visitor_delete(feature_stream->feature_visitor);
}

const GtNodeStreamClass* gt_feature_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtFeatureStream),
                                   feature_stream_free,
                                   feature_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_feature_stream_new(GtNodeStream *in_stream, GtFeatureIndex *fi)
{
  GtNodeStream *gs;
  GtFeatureStream *feature_stream;
  gs = gt_node_stream_create(gt_feature_stream_class(), false);
  feature_stream = feature_stream_cast(gs);
  feature_stream->in_stream = gt_node_stream_ref(in_stream);
  feature_stream->feature_visitor = gt_feature_visitor_new(fi);
  return gs;
}
