/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgt/transcript_counts.h>

struct TranscriptCounts {
  Array *exon_array_all,
        *exon_array_single,
        *exon_array_initial,
        *exon_array_internal,
        *exon_array_terminal;
};

TranscriptCounts* transcript_counts_new(Env *env)
{
  TranscriptCounts *tc = env_ma_calloc(env, 1, sizeof (TranscriptCounts));
  return tc;
}

Array* transcript_counts_get_all(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_all;
}

void transcript_counts_set_all(TranscriptCounts *tc, Array *counts)
{
  assert(tc && counts);
  tc->exon_array_all = counts;
}

Array* transcript_counts_get_single(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_single;
}

void transcript_counts_set_single(TranscriptCounts *tc, Array *counts)
{
  assert(tc && counts);
  tc->exon_array_single = counts;
}

Array* transcript_counts_get_initial(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_initial;
}

void transcript_counts_set_initial(TranscriptCounts *tc, Array *counts)
{
  assert(tc && counts);
  tc->exon_array_initial = counts;
}

Array* transcript_counts_get_internal(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_internal;
}

void transcript_counts_set_internal(TranscriptCounts *tc, Array *counts)
{
  assert(tc && counts);
  tc->exon_array_internal = counts;
}

Array* transcript_counts_get_terminal(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_terminal;
}

void transcript_counts_set_terminal(TranscriptCounts *tc, Array *counts)
{
  assert(tc && counts);
  tc->exon_array_terminal = counts;
}

void transcript_counts_delete(TranscriptCounts *tc, Env *env)
{
  if (!tc) return;
  array_delete(tc->exon_array_all, env);
  array_delete(tc->exon_array_single, env);
  array_delete(tc->exon_array_initial, env);
  array_delete(tc->exon_array_internal, env);
  array_delete(tc->exon_array_terminal, env);
  env_ma_free(tc, env);
}
