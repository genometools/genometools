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
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtext/genome_stream_rep.h"

GenomeStream* genome_stream_create(const GenomeStreamClass *gsc,
                                   bool ensure_sorting)
{
  GenomeStream *gs;
  assert(gsc && gsc->size);
  gs = ma_calloc(1, gsc->size);
  gs->c_class = gsc;
  gs->ensure_sorting = ensure_sorting;
  return gs;
}

GenomeStream* genome_stream_ref(GenomeStream *gs)
{
  assert(gs);
  gs->reference_count++;
  return gs;
}

void genome_stream_delete(GenomeStream *gs)
{
  if (!gs) return;
  if (gs->reference_count) {
    gs->reference_count--;
    return;
  }
  assert(gs->c_class);
  if (gs->c_class->free) gs->c_class->free(gs);
  genome_node_rec_delete(gs->buffer);
  ma_free(gs);
}

int genome_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Error *e)
{
  GenomeNode *new_node = NULL;
  int had_err = 0;
  assert(gs && gs->c_class && gs->c_class->next_tree);
  error_check(e);
  /* filling */
  if (!gs->buffer)
    had_err = gs->c_class->next_tree(gs, &gs->buffer, e);
  if (!had_err && gs->buffer)
    had_err = gs->c_class->next_tree(gs, &new_node, e);
#ifndef NDEBUG
  /* checking */
  if (!had_err && gs->ensure_sorting && gs->buffer && new_node) {
    assert(genome_node_compare(&gs->buffer, &new_node) <= 0);
  }
#endif
  /* serving */
  if (!had_err) {
    *gn = gs->buffer;
    gs->buffer = new_node;
  }
  return had_err;
}

bool genome_stream_is_sorted(GenomeStream *gs)
{
  assert(gs);
  return gs->ensure_sorting;
}

void* genome_stream_cast(UNUSED const GenomeStreamClass *gsc, GenomeStream *gs)
{
  assert(gsc && gs && gs->c_class == gsc);
  return gs;
}
