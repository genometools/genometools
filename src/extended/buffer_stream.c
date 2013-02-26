/*
  Copyright (c) 2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "extended/buffer_stream.h"
#include "extended/node_stream_api.h"

struct GtBufferStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtQueue *node_buffer;
  bool buffering;
};

#define buffer_stream_cast(NS)\
        gt_node_stream_cast(gt_buffer_stream_class(), NS)

static int buffer_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtBufferStream *bs;
  gt_error_check(err);
  bs = buffer_stream_cast(ns);
  if (bs->buffering) {
    int had_err = gt_node_stream_next(bs->in_stream, gn, err);
    if (!had_err && *gn)
      gt_queue_add(bs->node_buffer, gt_genome_node_ref(*gn));
    return had_err;
  }
  else {
    *gn = gt_queue_size(bs->node_buffer) ? gt_queue_get(bs->node_buffer) : NULL;
    return 0;
  }
}

static void buffer_stream_free(GtNodeStream *ns)
{
  GtBufferStream *bs = buffer_stream_cast(ns);
  while (gt_queue_size(bs->node_buffer))
    gt_genome_node_delete(gt_queue_get(bs->node_buffer));
  gt_queue_delete(bs->node_buffer);
  gt_node_stream_delete(bs->in_stream);
}

const GtNodeStreamClass* gt_buffer_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtBufferStream),
                                   buffer_stream_free,
                                   buffer_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_buffer_stream_new(GtNodeStream *in_stream)
{
  GtBufferStream *bs;
  GtNodeStream *ns;
  gt_assert(in_stream);
  ns = gt_node_stream_create(gt_buffer_stream_class(), false);
  bs = buffer_stream_cast(ns);
  bs->in_stream = gt_node_stream_ref(in_stream);
  bs->node_buffer = gt_queue_new();
  bs->buffering = true;
  return ns;
}

void gt_buffer_stream_dequeue(GtBufferStream *bs)
{
  gt_assert(bs);
  bs->buffering = false;
}
