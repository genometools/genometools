/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "transcript_bittabs.h"

struct TranscriptBittabs {
  Bittab *bittab_all,
         *bittab_single,
         *bittab_initial,
         *bittab_internal,
         *bittab_terminal;
};

TranscriptBittabs* transcript_bittabs_new(unsigned long size_all,
                                          unsigned long size_single,
                                          unsigned long size_initial,
                                          unsigned long size_internal,
                                          unsigned long size_terminal, Env *env)
{
  TranscriptBittabs *tc = env_ma_calloc(env, 1, sizeof (TranscriptBittabs));
  if (size_all) tc->bittab_all = bittab_new(size_all, env);
  if (size_single) tc->bittab_single = bittab_new(size_single, env);
  if (size_initial) tc->bittab_initial = bittab_new(size_initial, env);
  if (size_internal) tc->bittab_internal = bittab_new(size_internal, env);
  if (size_terminal) tc->bittab_terminal = bittab_new(size_all, env);
  return tc;
}

Bittab* transcript_bittabs_get_all(const TranscriptBittabs *tc)
{
  assert(tc);
  return tc->bittab_all;
}

Bittab* transcript_bittabs_get_single(const TranscriptBittabs *tc)
{
  assert(tc);
  return tc->bittab_single;
}

Bittab* transcript_bittabs_get_initial(const TranscriptBittabs *tc)
{
  assert(tc);
  return tc->bittab_initial;
}

Bittab* transcript_bittabs_get_internal(const TranscriptBittabs *tc)
{
  assert(tc);
  return tc->bittab_internal;
}

Bittab* transcript_bittabs_get_terminal(const TranscriptBittabs *tc)
{
  assert(tc);
  return tc->bittab_terminal;
}

void transcript_bittabs_delete(TranscriptBittabs *tc, Env *env)
{
  if (!tc) return;
  bittab_delete(tc->bittab_all, env);
  bittab_delete(tc->bittab_single, env);
  bittab_delete(tc->bittab_initial, env);
  bittab_delete(tc->bittab_internal, env);
  bittab_delete(tc->bittab_terminal, env);
  env_ma_free(tc, env);
}
