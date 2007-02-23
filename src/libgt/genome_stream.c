/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdarg.h>
#include "fptr.h"
#include "genome_stream_rep.h"
#include "xansi.h"

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

void genome_stream_delete(GenomeStream *gs, Env *env)
{
  if (!gs) return;
  assert(gs->c_class);
  if (gs->c_class->free) gs->c_class->free(gs, env);
  genome_node_delete(gs->last_node, env);
  env_ma_free(gs, env);
}

int genome_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  int has_err;
  assert(gs && gs->c_class && gs->c_class->next_tree);
  env_error_check(env);
  has_err = gs->c_class->next_tree(gs, gn, env);
  if (!has_err && *gn && gs->ensure_sorting) {
    assert(genome_node_tree_is_sorted(&gs->last_node, *gn, env));
  }
  return has_err;
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
