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

void genome_stream_class_init(GenomeStreamClass *gsc, size_t size, ...)
{
  va_list ap;
  Fptr func, meth, *mm;

  assert(gsc && size);

  gsc->size = size;
  va_start(ap, size);
  while ((func = va_arg(ap, Fptr))) {
    meth = va_arg(ap, Fptr);
    assert(meth);
    if (func == (Fptr) genome_stream_next_tree) {
      mm  = (Fptr*) &gsc->next_tree; *mm = meth;
    }
    else if (func == (Fptr) genome_stream_delete) {
      mm  = (Fptr*) &gsc->free; *mm = meth;
    }
    else assert(0);
  }
  va_end(ap);
}

GenomeStream* genome_stream_create(const GenomeStreamClass *gsc,
                                    bool ensure_sorting)
{
  GenomeStream *gs;
  assert(gsc && gsc->size);
  gs = xcalloc(1, gsc->size);
  gs->c_class = gsc;
  gs->ensure_sorting = ensure_sorting;
  return gs;
}

void genome_stream_delete(GenomeStream *gs)
{
  if (!gs) return;
  assert(gs->c_class);
  if (gs->c_class->free) gs->c_class->free(gs);
  genome_node_delete(gs->last_node);
  free(gs);
}

int genome_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Log *l,
                            Env *env)
{
  int has_err;
  assert(gs && gs->c_class && gs->c_class->next_tree);
  env_error_check(env);
  has_err = gs->c_class->next_tree(gs, gn, l, env);
  if (!has_err && *gn && gs->ensure_sorting) {
    assert(genome_node_tree_is_sorted(&gs->last_node, *gn));
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
