/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/ma.h"
#include "libgtext/transcript_evaluators.h"

struct TranscriptEvaluators {
  Evaluator *exon_evaluator_all,
            *exon_evaluator_single,
            *exon_evaluator_initial,
            *exon_evaluator_internal,
            *exon_evaluator_terminal;
};

TranscriptEvaluators* transcript_evaluators_new(void)
{
  TranscriptEvaluators *te = ma_malloc(sizeof (TranscriptEvaluators));
  te->exon_evaluator_all = evaluator_new();
  te->exon_evaluator_single = evaluator_new();
  te->exon_evaluator_initial = evaluator_new();
  te->exon_evaluator_internal = evaluator_new();
  te->exon_evaluator_terminal = evaluator_new();
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

void transcript_evaluators_delete(TranscriptEvaluators *te)
{
  if (!te) return;
  evaluator_delete(te->exon_evaluator_all);
  evaluator_delete(te->exon_evaluator_single);
  evaluator_delete(te->exon_evaluator_initial);
  evaluator_delete(te->exon_evaluator_internal);
  evaluator_delete(te->exon_evaluator_terminal);
  ma_free(te);
}
