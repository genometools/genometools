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
  TranscriptBittabs *tb = env_ma_calloc(env, 1, sizeof (TranscriptBittabs));
  if (size_all) tb->bittab_all = bittab_new(size_all, env);
  if (size_single) tb->bittab_single = bittab_new(size_single, env);
  if (size_initial) tb->bittab_initial = bittab_new(size_initial, env);
  if (size_internal) tb->bittab_internal = bittab_new(size_internal, env);
  if (size_terminal) tb->bittab_terminal = bittab_new(size_all, env);
  return tb;
}

Bittab* transcript_bittabs_get_all(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_all;
}

Bittab* transcript_bittabs_get_single(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_single;
}

Bittab* transcript_bittabs_get_initial(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_initial;
}

Bittab* transcript_bittabs_get_internal(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_internal;
}

Bittab* transcript_bittabs_get_terminal(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_terminal;
}

void transcript_bittabs_delete(TranscriptBittabs *tb, Env *env)
{
  if (!tb) return;
  bittab_delete(tb->bittab_all, env);
  bittab_delete(tb->bittab_single, env);
  bittab_delete(tb->bittab_initial, env);
  bittab_delete(tb->bittab_internal, env);
  bittab_delete(tb->bittab_terminal, env);
  env_ma_free(tb, env);
}
