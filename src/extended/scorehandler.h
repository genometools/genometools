/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef SCOREHANDLER_H
#define SCOREHANDLER_H
#include "core/alphabet.h"
#include "core/score_matrix.h"
#include "core/types_api.h"
#include "extended/alignment.h"

typedef struct GtScoreHandler GtScoreHandler;

/* A function to compare chars. */
typedef GtWord (Scorecomparefunc)(GtScoreHandler*, GtUchar, GtUchar);

GtScoreHandler* gt_scorehandler_new_DNA(GtWord matchscore,
                                        GtWord mismatchscore,
                                        GtWord gap_opening,
                                        GtWord gap_extension);

GtScoreHandler* gt_scorehandler_new_Protein(GtScoreMatrix *sm,
                                            GtWord gap_opening,
                                            GtWord gap_extension);

void gt_scorehandler_delete(GtScoreHandler *scorehandler);

GtWord gt_scorehandler_get_gap_opening(const GtScoreHandler *scorehandler);

GtWord gt_scorehandler_get_gapscore(const GtScoreHandler *scorehandler);

GtWord gt_scorehandler_get_replacement(GtScoreHandler *scorehandler,
                                       GtUchar a, GtUchar b);

GtScoreHandler *gt_scorehandler_get_costhandler(const GtScoreHandler
                                                                *scorehandler);

GtAlphabet *gt_scorehandler_get_alphabet(const GtScoreHandler *scorehandler);

void gt_scorehandler_change_score_to_cost(GtScoreHandler *scorehandler);

void gt_scorehandler_change_score_to_cost_without_costhandler(GtScoreHandler
                                                                 *scorehandler);

GtWord gt_scorehandler_eval_alignmentscore(const GtScoreHandler *scorehandler,
                                           const GtAlignment *alignment,
                                           const GtUchar *characters);

GtUchar* check_dna_sequence(const GtUchar *seq,
                            GtUword len, GtAlphabet *alphabet);
#endif
