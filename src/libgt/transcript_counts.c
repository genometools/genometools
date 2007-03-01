/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "transcript_counts.h"

struct TranscriptCounts {
  Array *exon_array_all,
        *exon_array_single,
        *exon_array_initial,
        *exon_array_internal,
        *exon_array_terminal;
};

TranscriptCounts* transcript_counts_new(Env *env)
{
  TranscriptCounts *te = env_ma_calloc(env, 1, sizeof (TranscriptCounts));
  return te;
}

Array* transcript_counts_get_all(const TranscriptCounts *te)
{
  assert(te);
  return te->exon_array_all;
}

void transcript_counts_set_all(TranscriptCounts *te, Array *counts)
{
  assert(te && counts);
  te->exon_array_all = counts;
}

Array* transcript_counts_get_single(const TranscriptCounts *te)
{
  assert(te);
  return te->exon_array_single;
}

void transcript_counts_set_single(TranscriptCounts *te, Array *counts)
{
  assert(te && counts);
  te->exon_array_single = counts;
}

Array* transcript_counts_get_initial(const TranscriptCounts *te)
{
  assert(te);
  return te->exon_array_initial;
}

void transcript_counts_set_initial(TranscriptCounts *te, Array *counts)
{
  assert(te && counts);
  te->exon_array_initial = counts;
}

Array* transcript_counts_get_internal(const TranscriptCounts *te)
{
  assert(te);
  return te->exon_array_internal;
}

void transcript_counts_set_internal(TranscriptCounts *te, Array *counts)
{
  assert(te && counts);
  te->exon_array_internal = counts;
}

Array* transcript_counts_get_terminal(const TranscriptCounts *te)
{
  assert(te);
  return te->exon_array_terminal;
}

void transcript_counts_set_terminal(TranscriptCounts *te, Array *counts)
{
  assert(te && counts);
  te->exon_array_terminal = counts;
}

void transcript_counts_delete(TranscriptCounts *te, Env *env)
{
  if (!te) return;
  array_delete(te->exon_array_all, env);
  array_delete(te->exon_array_single, env);
  array_delete(te->exon_array_initial, env);
  array_delete(te->exon_array_internal, env);
  array_delete(te->exon_array_terminal, env);
}
