/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/array.h"
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "extended/eof_node.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "extended/load_stream.h"

struct GtLoadStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  unsigned long idx;
  GtArray *nodes;
  bool full;
};

const GtNodeStreamClass* gt_load_stream_class(void);

#define gt_load_stream_cast(NS)\
        gt_node_stream_cast(gt_load_stream_class(), NS);

static int gt_load_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                               GtError *err)
{
  GtLoadStream *load_stream;
  GtGenomeNode *node, *eofn;
  int had_err = 0;
  gt_error_check(err);
  load_stream = gt_load_stream_cast(ns);

  if (!load_stream->full) {
    while (!(had_err = gt_node_stream_next(load_stream->in_stream, &node,
                                           err)) && node) {
      if ((eofn = gt_eof_node_try_cast(node)))
        gt_genome_node_delete(eofn); /* get rid of EOF nodes */
      else
        gt_array_add(load_stream->nodes, node);
    }
    if (!had_err) {
      load_stream->full = true;
    }
  }

  if (!had_err) {
    gt_assert(load_stream->full);
    if (load_stream->idx < gt_array_size(load_stream->nodes)) {
      *gn = *(GtGenomeNode**) gt_array_get(load_stream->nodes,
                                           load_stream->idx);
      load_stream->idx++;
      return 0;
    }
  }

  if (!had_err) {
    gt_array_reset(load_stream->nodes);
    *gn = NULL;
  }

  return had_err;
}

static void gt_load_stream_free(GtNodeStream *ns)
{
  unsigned long i;
  GtLoadStream *load_stream = gt_load_stream_cast(ns);
  for (i = load_stream->idx; i < gt_array_size(load_stream->nodes); i++) {
    gt_genome_node_delete(*(GtGenomeNode**)
                          gt_array_get(load_stream->nodes, i));
  }
  gt_array_delete(load_stream->nodes);
  gt_node_stream_delete(load_stream->in_stream);
}

const GtNodeStreamClass* gt_load_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLoadStream),
                                   gt_load_stream_free,
                                   gt_load_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_load_stream_new(GtNodeStream *in_stream)
{
  GtNodeStream *ns = gt_node_stream_create(gt_load_stream_class(), true);
  GtLoadStream *load_stream = gt_load_stream_cast(ns);
  gt_assert(in_stream);
  load_stream->in_stream = gt_node_stream_ref(in_stream);
  load_stream->full = false;
  load_stream->idx = 0;
  load_stream->nodes = gt_array_new(sizeof (GtGenomeNode*));
  return ns;
}
