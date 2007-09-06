/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/genome_stream_rep.h"

GenomeStream* genome_stream_create(const GenomeStreamClass *gsc,
                                   bool ensure_sorting, Env *env)
{
  GenomeStream *gs;
  assert(gsc && gsc->size);
  gs = env_ma_calloc(env, 1, gsc->size);
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

void genome_stream_delete(GenomeStream *gs, Env *env)
{
  if (!gs) return;
  if (gs->reference_count) {
    gs->reference_count--;
    return;
  }
  assert(gs->c_class);
  if (gs->c_class->free) gs->c_class->free(gs, env);
  genome_node_delete(gs->last_node, env);
  env_ma_free(gs, env);
}

int genome_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  int had_err;
  assert(gs && gs->c_class && gs->c_class->next_tree);
  env_error_check(env);
  had_err = gs->c_class->next_tree(gs, gn, env);
  if (!had_err && *gn && gs->ensure_sorting) {
    assert(genome_node_tree_is_sorted(&gs->last_node, *gn, env));
  }
  return had_err;
}

bool genome_stream_is_sorted(GenomeStream *gs)
{
  assert(gs);
  return gs->ensure_sorting;
}

void* genome_stream_cast(const GenomeStreamClass *gsc, GenomeStream *gs)
{
  assert(gsc && gs && gs->c_class == gsc);
  return gs;
}
