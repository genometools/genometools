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
#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/ma.h"
#include "core/score_matrix.h"

#include "extended/scoreHandler.h"
struct GtScoreHandler{
  GtWord           matchscore,
                   mismatchscore,
                   gap_opening,
                   gap_extension;
  GtScoreMatrix    *sm;
  GtAlphabet       *alphabet;
  Scorecomparefunc compare;
};

GtScoreHandler* gt_scorehandler_new_DNA(GtWord matchscore,
                                        GtWord mismatchscore,
                                        GtWord gap_opening,
                                        GtWord gap_extension)
{
  GtScoreHandler *scorehandler;
  scorehandler = gt_malloc(sizeof (*scorehandler));
  scorehandler->matchscore = matchscore;
  scorehandler->mismatchscore = mismatchscore;
  scorehandler->gap_opening = gap_opening;
  scorehandler->gap_extension = gap_extension;

  return scorehandler;
}

GtScoreHandler* gt_scorehandler_new_Protein(GtWord gap_opening,
                                            GtWord gap_extension,
                                            const char *path,
                                            GtError *err)
{
  GtScoreHandler *scorehandler;
  scorehandler = gt_malloc(sizeof (*scorehandler));
  scorehandler->sm = gt_score_matrix_new_read_protein(path,err);
  scorehandler->alphabet = gt_alphabet_new_protein();
  scorehandler->matchscore = GT_WORD_MAX;
  scorehandler->mismatchscore = GT_WORD_MAX;
  scorehandler->gap_opening = gap_opening;
  scorehandler->gap_extension = gap_extension;

  return scorehandler;
}

void gt_scorehandler_delete(GtScoreHandler *scorehandler)
{
  if (scorehandler != NULL)
  {
    gt_score_matrix_delete(scorehandler->sm);
    gt_alphabet_delete(scorehandler->alphabet);
    gt_free(scorehandler);
  }
}

GtWord gt_scorehandler_get_gap_opening(const GtScoreHandler *scorehandler)
{
  gt_assert(scorehandler != NULL);
  return scorehandler->gap_opening;
}

GtWord gt_scorehandler_get_gapscore(const GtScoreHandler *scorehandler)
{
  gt_assert(scorehandler != NULL);
  return scorehandler->gap_extension;
}

GtWord DNA_replacement(const GtScoreHandler *scorehandler, GtUchar a, GtUchar b)
{
  gt_assert(scorehandler != NULL && scorehandler->matchscore != GT_WORD_MAX &&
            scorehandler->matchscore != GT_WORD_MAX);
  /* matchscore or mismatchscore is only GT_WORD_MAX
   * if scorehandler is created with Proteins */

  if (ISSPECIAL(a) || ISSPECIAL(b) || a != b)
  {
    return scorehandler->mismatchscore;
  }

  return scorehandler->matchscore;
}

GtWord Protein_replacement(const GtScoreHandler *scorehandler,
                           GtUchar a, GtUchar b)
{
  GtUchar idx1, idx2;
  gt_assert(scorehandler != NULL);

//TODO: Wildcards?tolower?toupper?

  idx1 = gt_alphabet_encode(scorehandler->alphabet, a);
  idx2 = gt_alphabet_encode(scorehandler->alphabet, b);
  return gt_score_matrix_get_score(scorehandler->sm, idx1, idx2);
}
