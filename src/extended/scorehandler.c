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
#include <ctype.h>
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/ma.h"
#include "core/minmax.h"

#include "extended/scorehandler.h"
struct GtScoreHandler{
  GtWord           matchscore,
                   mismatchscore,
                   gap_opening,
                   gap_extension;
  GtScoreMatrix    *sm;
  GtAlphabet       *alphabet;
  Scorecomparefunc *compare;
  GtScoreHandler   *costhandler;/* only necessary for local-translations */
};

/* Scorecomparefunctions */
static GtWord DNA_replacement(const GtScoreHandler *scorehandler,
                              GtUchar a, GtUchar b);

static GtWord Protein_replacement(const GtScoreHandler *scorehandler,
                                  GtUchar a, GtUchar b);

GtScoreHandler* gt_scorehandler_new_DNA(GtWord matchscore,
                                        GtWord mismatchscore,
                                        GtWord gap_opening,
                                        GtWord gap_extension)
{
  GtScoreHandler *scorehandler;
  scorehandler = gt_malloc(sizeof (*scorehandler));
  scorehandler->sm = NULL;
  scorehandler->alphabet = gt_alphabet_new_dna();
  scorehandler->matchscore = matchscore;
  scorehandler->mismatchscore = mismatchscore;
  scorehandler->gap_opening = gap_opening;
  scorehandler->gap_extension = gap_extension;
  scorehandler->compare = (Scorecomparefunc*) DNA_replacement;
  scorehandler->costhandler = NULL;

  return scorehandler;
}

GtScoreHandler* gt_scorehandler_new_Protein(GtScoreMatrix *sm,
                                            GtWord gap_opening,
                                            GtWord gap_extension)
{
  GtScoreHandler *scorehandler;
  scorehandler = gt_malloc(sizeof (*scorehandler));
  scorehandler->sm = sm;
  scorehandler->alphabet = gt_alphabet_new_protein();
  scorehandler->matchscore = GT_WORD_MAX;
  scorehandler->mismatchscore = GT_WORD_MAX;
  scorehandler->gap_opening = gap_opening;
  scorehandler->gap_extension = gap_extension;
  scorehandler->compare = (Scorecomparefunc*) Protein_replacement;
  scorehandler->costhandler = NULL;

  return scorehandler;
}

void gt_scorehandler_delete(GtScoreHandler *scorehandler)
{
  if (scorehandler != NULL)
  {
    gt_score_matrix_delete(scorehandler->sm);
    gt_alphabet_delete(scorehandler->alphabet);
    gt_scorehandler_delete(scorehandler->costhandler);
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

static GtWord DNA_replacement(const GtScoreHandler *scorehandler,
                              GtUchar a, GtUchar b)
{
  gt_assert(scorehandler != NULL && scorehandler->matchscore != GT_WORD_MAX &&
            scorehandler->matchscore != GT_WORD_MAX);
  /* matchscore or mismatchscore is only GT_WORD_MAX
   * if scorehandler is created with Proteins */

  a = gt_alphabet_encode(scorehandler->alphabet, a);
  b = gt_alphabet_encode(scorehandler->alphabet, b);
  if (ISSPECIAL(a) || ISSPECIAL(b) || a != b)
  {
    return scorehandler->mismatchscore;
  }

  return scorehandler->matchscore;
}

static GtWord Protein_replacement(const GtScoreHandler *scorehandler,
                                  GtUchar a, GtUchar b)
{
  GtUchar idx1, idx2;
  gt_assert(scorehandler && scorehandler->sm && scorehandler->alphabet);

  idx1 = gt_alphabet_encode(scorehandler->alphabet, a);
  idx2 = gt_alphabet_encode(scorehandler->alphabet, b);

  int score=gt_score_matrix_get_score(scorehandler->sm, idx1, idx2);
  return score;
}

GtWord gt_scorehandler_get_replacement(GtScoreHandler *scorehandler,
                                       GtUchar a, GtUchar b)
{
  gt_assert(scorehandler != NULL);

  return scorehandler->compare(scorehandler, a, b);
}

GtScoreHandler *gt_scorehandler_get_costhandler(const GtScoreHandler
                                                *scorehandler)
{
  gt_assert(scorehandler != NULL);
  return scorehandler->costhandler;
}

GtAlphabet *gt_scorehandler_get_alphabet(const GtScoreHandler *scorehandler)
{
  gt_assert(scorehandler != NULL);
  return scorehandler->alphabet;
}

void gt_scorehandler_change_score_to_cost(GtScoreHandler *scorehandler)
{
  GtWord max, val, matchscore, mismatchscore, gap_extension, gap_opening;
  unsigned int i,j, dim;

  gt_assert(scorehandler);

  /* if costhandler exists already, translate scores to cost is not necessary */
  if (scorehandler->costhandler != NULL)
    return;

  if (scorehandler->sm == NULL)/* DNA */
  {
    max = MAX(MAX(GT_DIV2(scorehandler->matchscore+1),
                  GT_DIV2(scorehandler->mismatchscore+1)),
              MAX(1 + scorehandler->gap_extension,0));

    matchscore = 2 * max-scorehandler->matchscore;
    mismatchscore = 2 * max-scorehandler->mismatchscore;
    gap_extension = max - scorehandler->gap_extension;
    gap_opening = -scorehandler->gap_opening;

    scorehandler->costhandler = gt_scorehandler_new_DNA(matchscore,
                                                        mismatchscore,
                                                        gap_opening,
                                                        gap_extension);
  }
  else /* protein */
  {
    max = 0;
    dim = gt_score_matrix_get_dimension(scorehandler->sm);
    GtScoreMatrix *sm = gt_score_matrix_new(scorehandler->alphabet);

    for (i = 0; i < dim; i++)
    {
      for (j = 0; j < dim; j++)
      {
          val = gt_score_matrix_get_score(scorehandler->sm, i, j);
          if (val > max)
            max = val;
      }
    }

    max = MAX(GT_DIV2(max+1), 1 + scorehandler->gap_extension);
    for (i = 0; i < dim; i++)
    {
      for (j = 0; j < dim; j++)
      {
        /* translate */
        gt_score_matrix_set_score(sm, i, j, 2 * max-
                            gt_score_matrix_get_score(scorehandler->sm, i, j));
      }
    }
    gap_extension = max - scorehandler->gap_extension;
    gap_opening = -scorehandler->gap_opening;
    scorehandler->costhandler = gt_scorehandler_new_Protein(sm,
                                                            gap_opening,
                                                            gap_extension);
  }
}

void gt_scorehandler_change_score_to_cost_without_costhandler(GtScoreHandler
                                                                 *scorehandler)
{
  GtWord max, val;
  unsigned int i,j, dim;

  gt_assert(scorehandler);

  max = 0;
  dim = gt_score_matrix_get_dimension(scorehandler->sm);

  for (i = 0; i < dim; i++)
  {
    for (j = 0; j < dim; j++)
    {
        val = gt_score_matrix_get_score(scorehandler->sm, i, j);
        if (val > max)
          max = val;
    }
  }

  max = MAX(GT_DIV2(max+1), 1 + scorehandler->gap_extension);
  for (i = 0; i < dim; i++)
  {
    for (j = 0; j < dim; j++)
    {
      /* translate */
      gt_score_matrix_set_score(scorehandler->sm, i, j,  2* max-
                            gt_score_matrix_get_score(scorehandler->sm, i, j));
    }
  }
  scorehandler->gap_extension = max - scorehandler->gap_extension;
  scorehandler->gap_opening = -scorehandler->gap_opening;
}

GtWord gt_scorehandler_eval_alignmentscore(const GtScoreHandler *scorehandler,
                                           const GtAlignment *alignment,
                                           const GtUchar *characters)
{
  gt_assert(scorehandler && alignment && characters);

  if (gt_alphabet_is_protein(scorehandler->alphabet))
  {
    return gt_alignment_eval_with_affine_scorematrix(characters, alignment,
                                                   scorehandler->sm,
                                                   scorehandler->gap_opening,
                                                   scorehandler->gap_extension);
  }
  else
  {
    return gt_alignment_eval_with_mapped_affine_score(characters, alignment,
                                                   scorehandler->matchscore,
                                                   scorehandler->mismatchscore,
                                                   scorehandler->gap_opening,
                                                   scorehandler->gap_extension);
  }
}

GtUchar* check_dna_sequence(const GtUchar *seq, GtUword len,
                            GtAlphabet *alphabet)
{
  GtUword i;
  GtUchar *low_seq;

  low_seq = gt_malloc(sizeof(*low_seq)*(len+1));
  for (i = 0; i < len; i++)
  {
    low_seq[i] = tolower((int)seq[i]);
    if (!gt_alphabet_valid_input(alphabet, seq[i]) ||
       ISSPECIAL(gt_alphabet_encode(alphabet, seq[i])))
    {
      /* gt_assert(false);*/
      gt_free(low_seq);
      return NULL;
    }
  }
  low_seq[i] = '\0';

  return low_seq;
}
