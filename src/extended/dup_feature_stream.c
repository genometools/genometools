/*
  Copyright (c) 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "extended/dup_feature_stream.h"
#include "extended/dup_feature_visitor.h"

struct GtDupFeatureStream{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *dup_feature_visitor;
};

#define gt_dup_feature_stream_cast(GS)\
        gt_node_stream_cast(gt_dup_feature_stream_class(), GS)

static int dup_feature_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                   GtError *err)
{
  GtDupFeatureStream *ais;
  int had_err;
  gt_error_check(err);
  ais = gt_dup_feature_stream_cast(ns);
  had_err = gt_node_stream_next(ais->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, ais->dup_feature_visitor, err);
  return had_err;
}

static void dup_feature_stream_free(GtNodeStream *ns)
{
  GtDupFeatureStream *ais = gt_dup_feature_stream_cast(ns);
  gt_node_visitor_delete(ais->dup_feature_visitor);
  gt_node_stream_delete(ais->in_stream);
}

const GtNodeStreamClass* gt_dup_feature_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
   nsc = gt_node_stream_class_new(sizeof (GtDupFeatureStream),
                                  dup_feature_stream_free,
                                  dup_feature_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_dup_feature_stream_new(GtNodeStream *in_stream,
                                        const char *dest_type,
                                        const char *source_type)
{
  GtNodeStream *ns = gt_node_stream_create(gt_dup_feature_stream_class(),
                                           false);
  GtDupFeatureStream *ais = gt_dup_feature_stream_cast(ns);
  gt_assert(in_stream);
  ais->in_stream = gt_node_stream_ref(in_stream);
  ais->dup_feature_visitor = gt_dup_feature_visitor_new(dest_type, source_type);
  return ns;
}
