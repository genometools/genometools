/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/assert_api.h"
#include "extended/node_stream_rep.h"
#include "extended/inter_feature_stream.h"
#include "extended/inter_feature_visitor.h"

struct GtInterFeatureStream{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *inter_feature_visitor;
};

#define gt_inter_feature_stream_cast(GS)\
        gt_node_stream_cast(gt_inter_feature_stream_class(), GS)

static int inter_feature_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                        GtError *err)
{
  GtInterFeatureStream *ais;
  int had_err;
  gt_error_check(err);
  ais = gt_inter_feature_stream_cast(ns);
  had_err = gt_node_stream_next(ais->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, ais->inter_feature_visitor, err);
  return had_err;
}

static void inter_feature_stream_free(GtNodeStream *ns)
{
  GtInterFeatureStream *ais = gt_inter_feature_stream_cast(ns);
  gt_node_visitor_delete(ais->inter_feature_visitor);
  gt_node_stream_delete(ais->in_stream);
}

const GtNodeStreamClass* gt_inter_feature_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
   nsc = gt_node_stream_class_new(sizeof (GtInterFeatureStream),
                                  inter_feature_stream_free,
                                  inter_feature_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_inter_feature_stream_new(GtNodeStream *in_stream,
                                             const char *outside_type,
                                             const char *inter_type)
{
  GtNodeStream *ns = gt_node_stream_create(gt_inter_feature_stream_class(),
                                           false);
  GtInterFeatureStream *ais = gt_inter_feature_stream_cast(ns);
  gt_assert(in_stream);
  ais->in_stream = gt_node_stream_ref(in_stream);
  ais->inter_feature_visitor = gt_inter_feature_visitor_new(outside_type,
                                                            inter_type);
  return ns;
}
