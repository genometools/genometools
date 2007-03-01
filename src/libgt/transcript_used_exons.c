/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "range.h"
#include "transcript_used_exons.h"

struct TranscriptUsedExons {
  Dlist *used_exons_all,
        *used_exons_single,
        *used_exons_initial,
        *used_exons_internal,
        *used_exons_terminal;
};

TranscriptUsedExons* transcript_used_exons_new(Env *env)
{
  TranscriptUsedExons *tue = env_ma_malloc(env, sizeof (TranscriptUsedExons));
  tue->used_exons_all = dlist_new((Compare) range_compare_ptr, env);
  tue->used_exons_single = dlist_new((Compare) range_compare_ptr, env);
  tue->used_exons_initial = dlist_new((Compare) range_compare_ptr, env);
  tue->used_exons_internal = dlist_new((Compare) range_compare_ptr, env);
  tue->used_exons_terminal = dlist_new((Compare) range_compare_ptr, env);
  return tue;
}

Dlist* transcript_used_exons_get_all(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_all;
}

Dlist* transcript_used_exons_get_single(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_single;
}

Dlist* transcript_used_exons_get_initial(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_initial;
}

Dlist* transcript_used_exons_get_internal(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_internal;
}

Dlist* transcript_used_exons_get_terminal(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_terminal;
}

static void used_dlist_delete(Dlist *used_list, Env *env)
{
  Dlistelem *dlistelem;
  for (dlistelem = dlist_first(used_list); dlistelem != NULL;
       dlistelem = dlistelem_next(dlistelem)) {
    env_ma_free(dlistelem_get_data(dlistelem), env);
  }
  dlist_delete(used_list, env);
}

void transcript_used_exons_delete(TranscriptUsedExons *tue, Env *env)
{
  if (!tue) return;
  used_dlist_delete(tue->used_exons_all, env);
  used_dlist_delete(tue->used_exons_single, env);
  used_dlist_delete(tue->used_exons_initial, env);
  used_dlist_delete(tue->used_exons_internal, env);
  used_dlist_delete(tue->used_exons_terminal, env);
  env_ma_free(tue, env);
}
