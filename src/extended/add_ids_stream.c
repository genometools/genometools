/*
  Copyright (c) 2010, 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "extended/add_ids_stream.h"
#include "extended/add_ids_visitor.h"
#include "extended/feature_node.h"
#include "extended/node_stream_api.h"

struct GtAddIDsStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *add_ids_visitor; /* the actual work is done in the visitor */
  bool add_ids;
};

#define gt_add_ids_stream_cast(GS)\
        gt_node_stream_cast(gt_add_ids_stream_class(), GS);

static int add_ids_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                               GtError *err)
{
  GtAddIDsStream *ais;
  int had_err;
  gt_error_check(err);
  ais = gt_add_ids_stream_cast(ns);

  /* we still have nodes in the buffer */
  if (gt_add_ids_visitor_node_buffer_size(ais->add_ids_visitor)) {
    /* return one of them */
    *gn = gt_add_ids_visitor_get_node(ais->add_ids_visitor);
    return 0;
  }

  /* handle disabled stream */
  if (!ais->add_ids) {
    had_err = gt_node_stream_next(ais->in_stream, gn, err);
    if (had_err) {
      /* we own the node -> delete it */
      gt_genome_node_delete(*gn);
      *gn = NULL;
    }
    return had_err;
  }

  /* no nodes in the buffer -> get new nodes */
  while (!(had_err = gt_node_stream_next(ais->in_stream, gn, err)) && *gn) {
    gt_assert(*gn && !had_err);
    had_err = gt_genome_node_accept(*gn, ais->add_ids_visitor, err);
    if (had_err) {
      /* we own the node -> delete it */
      gt_genome_node_delete(*gn);
      *gn = NULL;
      break;
    }
    if (gt_add_ids_visitor_node_buffer_size(ais->add_ids_visitor)) {
      *gn = gt_add_ids_visitor_get_node(ais->add_ids_visitor);
      return 0;
    }
  }

  if (!had_err) {
    gt_add_ids_visitor_finalize(ais->add_ids_visitor);
    if (gt_add_ids_visitor_node_buffer_size(ais->add_ids_visitor)) {
      *gn = gt_add_ids_visitor_get_node(ais->add_ids_visitor);
      return 0;
    }
  }

  /* either we have an error or no new node */
  gt_assert(had_err || !*gn);
  return had_err;
}

static void add_ids_stream_free(GtNodeStream *ns)
{
  GtAddIDsStream *ais = gt_add_ids_stream_cast(ns);
  gt_node_visitor_delete(ais->add_ids_visitor);
  gt_node_stream_delete(ais->in_stream);
}

const GtNodeStreamClass* gt_add_ids_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtAddIDsStream),
                                   add_ids_stream_free,
                                   add_ids_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_add_ids_stream_new(GtNodeStream *in_stream)
{
  GtAddIDsStream *add_ids_stream;
  GtNodeStream *ns;
  gt_assert(in_stream);
  ns = gt_node_stream_create(gt_add_ids_stream_class(), false);
  add_ids_stream = gt_add_ids_stream_cast(ns);
  add_ids_stream->in_stream = gt_node_stream_ref(in_stream);
  add_ids_stream->add_ids_visitor =
    gt_add_ids_visitor_new(gt_node_stream_is_sorted(in_stream));
  add_ids_stream->add_ids = true;
  return ns;
}

void gt_add_ids_stream_disable(GtNodeStream *ns)
{
  GtAddIDsStream *add_ids_stream = gt_add_ids_stream_cast(ns);
  add_ids_stream->add_ids = false;
}
