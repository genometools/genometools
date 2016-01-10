/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
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
#include <ctype.h>
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "extended/scorehandler.h"

struct GtScoreHandler {
  GtWord matchscore,
         mismatchscore,
         gap_opening,
         gap_extension;
  GtScoreMatrix *scorematrix;
  bool mappedsequence,
       downcase;
};

/* Scorecomparefunctions */

GtScoreHandler* gt_scorehandler_new(GtWord matchscore,
                                    GtWord mismatchscore,
                                    GtWord gap_opening,
                                    GtWord gap_extension)
{
  GtScoreHandler *scorehandler = gt_malloc(sizeof *scorehandler);

  scorehandler->scorematrix = NULL;
  scorehandler->mappedsequence = true;
  scorehandler->matchscore = matchscore;
  scorehandler->mismatchscore = mismatchscore;
  scorehandler->gap_opening = gap_opening;
  scorehandler->gap_extension = gap_extension;
  scorehandler->downcase = false;
  return scorehandler;
}

void gt_scorehandler_add_scorematrix(GtScoreHandler *scorehandler,
                                     GtScoreMatrix *scorematrix)
{
  gt_assert(scorehandler != NULL);
  scorehandler->scorematrix = scorematrix;
}

void gt_scorehandler_plain(GtScoreHandler *scorehandler)
{
  gt_assert(scorehandler != NULL);
  scorehandler->mappedsequence = false;
}

void gt_scorehandler_downcase(GtScoreHandler *scorehandler)
{
  gt_assert(scorehandler != NULL);
  scorehandler->downcase = true;
}

void gt_scorehandler_delete(GtScoreHandler *scorehandler)
{
  if (scorehandler != NULL)
  {
    gt_score_matrix_delete(scorehandler->scorematrix);
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

GtWord gt_scorehandler_get_replacement(const GtScoreHandler *scorehandler,
                                       GtUchar a, GtUchar b)
{
  gt_assert(scorehandler != NULL);
  if (scorehandler->scorematrix == NULL)
  {
    if (scorehandler->mappedsequence)
    {
      return ISSPECIAL(a) || ISSPECIAL(b) || a != b
               ? scorehandler->mismatchscore
               : scorehandler->matchscore;
    }
    if (scorehandler->downcase)
    {
      a = tolower((int) a);
      b = tolower((int) b);
    }
    return a != b ? scorehandler->mismatchscore
                  : scorehandler->matchscore;
  }
  gt_assert(scorehandler->mappedsequence);
  return gt_score_matrix_get_score(scorehandler->scorematrix,a,b);
}

GtScoreHandler *gt_scorehandler2costhandler(const GtScoreHandler *scorehandler)
{
  GtScoreHandler *costhandler;

  gt_assert(scorehandler != NULL);
  if (scorehandler->scorematrix == NULL)
  {
    GtWord matchscore, mismatchscore, gap_extension, gap_opening,
           maxscore = MAX(MAX(GT_DIV2(scorehandler->matchscore+1),
                         GT_DIV2(scorehandler->mismatchscore+1)),
                     MAX(1 + scorehandler->gap_extension,0));

    matchscore = 2 * maxscore - scorehandler->matchscore;
    mismatchscore = 2 * maxscore - scorehandler->mismatchscore;
    gap_extension = maxscore - scorehandler->gap_extension;
    gap_opening = -scorehandler->gap_opening;
    costhandler = gt_scorehandler_new(matchscore,
                                      mismatchscore,
                                      gap_opening,
                                      gap_extension);
    if (!scorehandler->mappedsequence)
    {
      gt_scorehandler_plain(costhandler);
    }
  } else
  {
    int maxscore;
    GtWord gap_extension, gap_opening;
    unsigned int i, j,
                 dim = gt_score_matrix_get_dimension(scorehandler->scorematrix);
    GtScoreMatrix *costmatrix
      = gt_score_matrix_clone_empty(scorehandler->scorematrix);

    for (maxscore = 0, i = 0; i < dim; i++)
    {
      for (j = 0; j < dim; j++)
      {
        int val = gt_score_matrix_get_score(scorehandler->scorematrix, i, j);

        if (val > maxscore)
        {
          maxscore = val;
        }
      }
    }
    maxscore = MAX(GT_DIV2(maxscore+1), 1 + scorehandler->gap_extension);
    for (i = 0; i < dim; i++)
    {
      for (j = 0; j < dim; j++)
      {
        /* translate */
        int score = gt_score_matrix_get_score(scorehandler->scorematrix,i,j);
        gt_score_matrix_set_score(costmatrix, i, j, 2 * maxscore - score);
      }
    }
    gap_extension = maxscore - scorehandler->gap_extension;
    gap_opening = -scorehandler->gap_opening;
    costhandler = gt_scorehandler_new( 0,0, gap_opening, gap_extension);
    gt_scorehandler_add_scorematrix(costhandler,costmatrix);
  }
  return costhandler;
}

GtWord gt_scorehandler_eval_alignmentscore(const GtScoreHandler *scorehandler,
                                           const GtAlignment *alignment,
                                           const GtUchar *characters)
{
  gt_assert(scorehandler && alignment && characters);
  if (scorehandler->scorematrix != NULL)
  {
    return gt_alignment_eval_with_affine_scorematrix(characters, alignment,
                                                   scorehandler->scorematrix,
                                                   scorehandler->gap_opening,
                                                   scorehandler->gap_extension);
  } else
  {
    return gt_alignment_eval_with_mapped_affine_score(characters, alignment,
                                                   scorehandler->matchscore,
                                                   scorehandler->mismatchscore,
                                                   scorehandler->gap_opening,
                                                   scorehandler->gap_extension);
  }
}
