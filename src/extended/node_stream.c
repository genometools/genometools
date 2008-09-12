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

#include <assert.h>
#include <stdarg.h>
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/node_stream_rep.h"

GtNodeStream* gt_node_stream_create(const GtNodeStreamClass *gsc,
                                   bool ensure_sorting)
{
  GtNodeStream *gs;
  assert(gsc && gsc->size);
  gs = gt_calloc(1, gsc->size);
  gs->c_class = gsc;
  gs->ensure_sorting = ensure_sorting;
  return gs;
}

GtNodeStream* gt_node_stream_ref(GtNodeStream *gs)
{
  assert(gs);
  gs->reference_count++;
  return gs;
}

void gt_node_stream_delete(GtNodeStream *gs)
{
  if (!gs) return;
  if (gs->reference_count) {
    gs->reference_count--;
    return;
  }
  assert(gs->c_class);
  if (gs->c_class->free) gs->c_class->free(gs);
  gt_genome_node_rec_delete(gs->buffer);
  gt_free(gs);
}

int gt_node_stream_next(GtNodeStream *gs, GtGenomeNode **gn, GtError *err)
{
  GtGenomeNode *new_node = NULL;
  int had_err = 0;
  assert(gs && gs->c_class && gs->c_class->next_tree);
  gt_error_check(err);
  /* filling */
  if (!gs->buffer)
    had_err = gs->c_class->next_tree(gs, &gs->buffer, err);
  if (!had_err && gs->buffer)
    had_err = gs->c_class->next_tree(gs, &new_node, err);
#ifndef NDEBUG
  /* checking */
  if (!had_err && gs->ensure_sorting && gs->buffer && new_node) {
    assert(gt_genome_node_compare(&gs->buffer, &new_node) <= 0);
  }
#endif
  /* serving */
  if (!had_err) {
    *gn = gs->buffer;
    gs->buffer = new_node;
  }
  return had_err;
}

bool gt_node_stream_is_sorted(GtNodeStream *gs)
{
  assert(gs);
  return gs->ensure_sorting;
}

void* gt_node_stream_cast(GT_UNUSED const GtNodeStreamClass *gsc,
                         GtNodeStream *gs)
{
  assert(gsc && gs && gs->c_class == gsc);
  return gs;
}
