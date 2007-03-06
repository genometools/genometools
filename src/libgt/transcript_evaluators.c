/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgt/transcript_evaluators.h>

struct TranscriptEvaluators {
  Evaluator *exon_evaluator_all,
            *exon_evaluator_single,
            *exon_evaluator_initial,
            *exon_evaluator_internal,
            *exon_evaluator_terminal;
};

TranscriptEvaluators* transcript_evaluators_new(Env *env)
{
  TranscriptEvaluators *te = env_ma_malloc(env, sizeof (TranscriptEvaluators));
  te->exon_evaluator_all = evaluator_new(env);
  te->exon_evaluator_single = evaluator_new(env);
  te->exon_evaluator_initial = evaluator_new(env);
  te->exon_evaluator_internal = evaluator_new(env);
  te->exon_evaluator_terminal = evaluator_new(env);
  return te;
}

Evaluator* transcript_evaluators_get_all(const TranscriptEvaluators *te)
{
  assert(te);
  return te->exon_evaluator_all;
}

Evaluator* transcript_evaluators_get_single(const TranscriptEvaluators *te)
{
  assert(te);
  return te->exon_evaluator_single;
}

Evaluator* transcript_evaluators_get_initial(const TranscriptEvaluators *te)
{
  assert(te);
  return te->exon_evaluator_initial;
}

Evaluator* transcript_evaluators_get_internal(const TranscriptEvaluators *te)
{
  assert(te);
  return te->exon_evaluator_internal;
}

Evaluator* transcript_evaluators_get_terminal(const TranscriptEvaluators *te)
{
  assert(te);
  return te->exon_evaluator_terminal;
}

void transcript_evaluators_add_actuals(const TranscriptEvaluators *evaluators,
                                       const TranscriptExons *exons)
{
  assert(evaluators && exons);
  evaluator_add_actual(evaluators->exon_evaluator_all,
                       array_size(transcript_exons_get_all(exons)));
  evaluator_add_actual(evaluators->exon_evaluator_single,
                       array_size(transcript_exons_get_single(exons)));
  evaluator_add_actual(evaluators->exon_evaluator_initial,
                       array_size(transcript_exons_get_initial(exons)));
  evaluator_add_actual(evaluators->exon_evaluator_internal,
                       array_size(transcript_exons_get_internal(exons)));
  evaluator_add_actual(evaluators->exon_evaluator_terminal,
                       array_size(transcript_exons_get_terminal(exons)));
}

void transcript_evaluators_delete(TranscriptEvaluators *te, Env *env)
{
  if (!te) return;
  evaluator_delete(te->exon_evaluator_all, env);
  evaluator_delete(te->exon_evaluator_single, env);
  evaluator_delete(te->exon_evaluator_initial, env);
  evaluator_delete(te->exon_evaluator_internal, env);
  evaluator_delete(te->exon_evaluator_terminal, env);
  env_ma_free(te, env);
}
