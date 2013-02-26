/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "extended/feature_node.h"
#include "extended/tidy_region_node_stream.h"
#include "extended/tidy_region_node_visitor.h"
#include "extended/node_stream_api.h"

struct GtTidyRegionNodeStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *tidy_region_node_visitor;
  bool queued;
};

#define gt_tidy_region_node_stream_cast(GS)\
        gt_node_stream_cast(gt_tidy_region_node_stream_class(), GS);

static int tidy_region_node_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                               GtError *err)
{
  GtTidyRegionNodeStream *ais;
  int had_err = 0;
  gt_error_check(err);
  ais = gt_tidy_region_node_stream_cast(ns);

  if (!ais->queued) {
    /* no nodes in the buffer -> get new nodes */
    while (!(had_err = gt_node_stream_next(ais->in_stream, gn, err)) && *gn) {
      gt_assert(*gn && !had_err);
      had_err = gt_genome_node_accept(*gn, ais->tidy_region_node_visitor, err);
      if (had_err) {
        /* we own the node -> delete it */
        gt_genome_node_delete(*gn);
        *gn = NULL;
        break;
      }
    }
    if (!had_err) {
      ais->queued = true;
    }
  }
  if (!had_err) {
    if (gt_tidy_region_node_visitor_node_buffer_size(
                                               ais->tidy_region_node_visitor)) {
      *gn = gt_tidy_region_node_visitor_get_node(ais->tidy_region_node_visitor);
    }
  }

  return had_err;
}

static void tidy_region_node_stream_free(GtNodeStream *ns)
{
  GtTidyRegionNodeStream *ais = gt_tidy_region_node_stream_cast(ns);
  gt_node_visitor_delete(ais->tidy_region_node_visitor);
  gt_node_stream_delete(ais->in_stream);
}

const GtNodeStreamClass* gt_tidy_region_node_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtTidyRegionNodeStream),
                                   tidy_region_node_stream_free,
                                   tidy_region_node_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_tidy_region_node_stream_new(GtNodeStream *in_stream)
{
  GtTidyRegionNodeStream *tidy_region_node_stream;
  GtNodeStream *ns;
  gt_assert(in_stream);
  ns = gt_node_stream_create(gt_tidy_region_node_stream_class(), false);
  tidy_region_node_stream = gt_tidy_region_node_stream_cast(ns);
  tidy_region_node_stream->in_stream = gt_node_stream_ref(in_stream);
  tidy_region_node_stream->tidy_region_node_visitor =
                                              gt_tidy_region_node_visitor_new();
  tidy_region_node_stream->queued = false;
  return ns;
}
