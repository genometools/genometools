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
#include "extended/cds_check_stream.h"
#include "extended/cds_check_visitor.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"

struct GtCDSCheckStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *cds_check_visitor;
};

#define cds_check_stream_cast(NS)\
        gt_node_stream_cast(gt_cds_check_stream_class(), NS)

static int cds_check_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                 GtError *err)
{
  GtCDSCheckStream *cs = cds_check_stream_cast(ns);
  int had_err;
  gt_error_check(err);
  cs = cds_check_stream_cast(ns);
  had_err = gt_node_stream_next(cs->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, cs->cds_check_visitor, err);
  if (had_err) {
    /* we own the node -> delete it */
    gt_genome_node_delete(*gn);
    *gn = NULL;
  }
  return had_err;
}

static void cds_check_stream_free(GtNodeStream *ns)
{
  GtCDSCheckStream *cs = cds_check_stream_cast(ns);
  gt_node_stream_delete(cs->in_stream);
  gt_node_visitor_delete(cs->cds_check_visitor);
}

const GtNodeStreamClass* gt_cds_check_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtCDSCheckStream),
                                   cds_check_stream_free,
                                   cds_check_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_cds_check_stream_new(GtNodeStream *in_stream)
{
  GtNodeStream *ns = gt_node_stream_create(gt_cds_check_stream_class(), false);
  GtCDSCheckStream *cs = cds_check_stream_cast(ns);
  cs->in_stream = gt_node_stream_ref(in_stream);
  cs->cds_check_visitor = gt_cds_check_visitor_new();
  return ns;
}

void gt_cds_check_stream_enable_tidy_mode(GtCDSCheckStream *cs)
{
  gt_assert(cs);
  gt_cds_check_visitor_enable_tidy_mode((GtCDSCheckVisitor*)
                                        cs->cds_check_visitor);
}
