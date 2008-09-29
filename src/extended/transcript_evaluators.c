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

#include "core/ma.h"
#include "extended/transcript_evaluators.h"

struct GtTranscriptEvaluators {
  GtEvaluator *exon_evaluator_all,
            *exon_evaluator_single,
            *exon_evaluator_initial,
            *exon_evaluator_internal,
            *exon_evaluator_terminal;
};

GtTranscriptEvaluators* gt_transcript_evaluators_new(void)
{
  GtTranscriptEvaluators *te = gt_malloc(sizeof (GtTranscriptEvaluators));
  te->exon_evaluator_all = gt_evaluator_new();
  te->exon_evaluator_single = gt_evaluator_new();
  te->exon_evaluator_initial = gt_evaluator_new();
  te->exon_evaluator_internal = gt_evaluator_new();
  te->exon_evaluator_terminal = gt_evaluator_new();
  return te;
}

GtEvaluator* gt_transcript_evaluators_get_all(const GtTranscriptEvaluators *te)
{
  gt_assert(te);
  return te->exon_evaluator_all;
}

GtEvaluator* gt_transcript_evaluators_get_single(const GtTranscriptEvaluators
                                                 *te)
{
  gt_assert(te);
  return te->exon_evaluator_single;
}

GtEvaluator* gt_transcript_evaluators_get_initial(const GtTranscriptEvaluators
                                                  *te)
{
  gt_assert(te);
  return te->exon_evaluator_initial;
}

GtEvaluator* gt_transcript_evaluators_get_internal(const GtTranscriptEvaluators
                                                   *te)
{
  gt_assert(te);
  return te->exon_evaluator_internal;
}

GtEvaluator* gt_transcript_evaluators_get_terminal(const GtTranscriptEvaluators
                                                   *te)
{
  gt_assert(te);
  return te->exon_evaluator_terminal;
}

void gt_transcript_evaluators_add_actuals(const GtTranscriptEvaluators
                                                                *evaluators,
                                          const GtTranscriptExons *exons)
{
  gt_assert(evaluators && exons);
  gt_evaluator_add_actual(evaluators->exon_evaluator_all,
                       gt_array_size(gt_transcript_exons_get_all(exons)));
  gt_evaluator_add_actual(evaluators->exon_evaluator_single,
                       gt_array_size(gt_transcript_exons_get_single(exons)));
  gt_evaluator_add_actual(evaluators->exon_evaluator_initial,
                       gt_array_size(gt_transcript_exons_get_initial(exons)));
  gt_evaluator_add_actual(evaluators->exon_evaluator_internal,
                       gt_array_size(gt_transcript_exons_get_internal(exons)));
  gt_evaluator_add_actual(evaluators->exon_evaluator_terminal,
                       gt_array_size(gt_transcript_exons_get_terminal(exons)));
}

void gt_transcript_evaluators_delete(GtTranscriptEvaluators *te)
{
  if (!te) return;
  gt_evaluator_delete(te->exon_evaluator_all);
  gt_evaluator_delete(te->exon_evaluator_single);
  gt_evaluator_delete(te->exon_evaluator_initial);
  gt_evaluator_delete(te->exon_evaluator_internal);
  gt_evaluator_delete(te->exon_evaluator_terminal);
  gt_free(te);
}
