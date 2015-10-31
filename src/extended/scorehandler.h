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
#include "core/score_matrix.h"
#include "core/types_api.h"
#include "extended/alignment.h"

typedef struct GtScoreHandler GtScoreHandler;

GtScoreHandler* gt_scorehandler_new(GtWord matchscore,
                                    GtWord mismatchscore,
                                    GtWord gap_opening,
                                    GtWord gap_extension);

void gt_scorehandler_add_scorematrix(GtScoreHandler *scorehandler,
                                     GtScoreMatrix *scorematrix);

void gt_scorehandler_plain(GtScoreHandler *scorehandler);

void gt_scorehandler_downcase(GtScoreHandler *scorehandler);

void gt_scorehandler_delete(GtScoreHandler *scorehandler);

GtWord gt_scorehandler_get_gap_opening(const GtScoreHandler *scorehandler);

GtWord gt_scorehandler_get_gapscore(const GtScoreHandler *scorehandler);

GtWord gt_scorehandler_get_replacement(const GtScoreHandler *scorehandler,
                                       GtUchar a, GtUchar b);

GtScoreHandler *gt_scorehandler2costhandler(const GtScoreHandler *scorehandler);

GtWord gt_scorehandler_eval_alignmentscore(const GtScoreHandler *scorehandler,
                                           const GtAlignment *alignment,
                                           const GtUchar *characters);

#endif
