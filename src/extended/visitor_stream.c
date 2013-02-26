/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "extended/feature_node.h"
#include "extended/node_stream_api.h"
#include "extended/visitor_stream.h"

struct GtVisitorStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *visitor;
};

#define visitor_stream_cast(GS)\
        gt_node_stream_cast(gt_visitor_stream_class(), GS)

static int visitor_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                               GtError *err)
{
  GtVisitorStream *visitor_stream;
  int had_err;
  gt_error_check(err);
  visitor_stream = visitor_stream_cast(ns);
  had_err = gt_node_stream_next(visitor_stream->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, visitor_stream->visitor, err);
  if (had_err) {
    /* we own the node -> delete it */
    gt_genome_node_delete(*gn);
    *gn = NULL;
  }
  return had_err;
}

static void visitor_stream_free(GtNodeStream *ns)
{
  GtVisitorStream *visitor_stream = visitor_stream_cast(ns);
  gt_node_visitor_delete(visitor_stream->visitor);
  gt_node_stream_delete(visitor_stream->in_stream);
}

const GtNodeStreamClass* gt_visitor_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtVisitorStream),
                                   visitor_stream_free,
                                   visitor_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_visitor_stream_new(GtNodeStream *in_stream,
                                    GtNodeVisitor *visitor)
{
  GtVisitorStream *visitor_stream;
  GtNodeStream *ns;
  ns = gt_node_stream_create(gt_visitor_stream_class(),
                             gt_node_stream_is_sorted(in_stream));
  visitor_stream = visitor_stream_cast(ns);
  visitor_stream->in_stream = gt_node_stream_ref(in_stream);
  visitor_stream->visitor = visitor;
  return ns;
}
