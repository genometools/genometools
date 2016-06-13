/*
  Copyright (c) 2016 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2016 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include "core/alphabet.h"
#include "core/ma.h"
#include "core/types_api.h"
#include "match/karlin_altschul_stat.h"
#include "extended/scorehandler.h"

/* TODO:reference, analog to blast */

struct GtKarlinAltschulStat
{
  double lambda,
         K,
         logK,
         H,
         alpha_div_lambda, /* TODO: evt alpha und beta nicht hier speichern*/
         beta;
};

typedef struct{
  double *sprob;
  GtWord low_score,
         high_score;
         /*min_score,
         max_score;*/
} ScoringFrequency;

typedef struct{
  char   ch;
  double p;
} LetterProb;

/* provisional solution */
static LetterProb nt_prob[] = {
  { 'A', 25.00 },
  { 'C', 25.00 },
  { 'G', 25.00 },
  { 'T', 25.00 }
};

GtKarlinAltschulStat *gt_ka_new(void)
{
  GtKarlinAltschulStat *ka;
  ka = gt_malloc(sizeof (GtKarlinAltschulStat));
  ka->lambda = 0;
  ka->K = 0;
  ka->logK = 0;
  ka->H = 0;
  ka->alpha_div_lambda = 0;
  ka->beta = 0;
  return ka;
}

void gt_ka_delete(GtKarlinAltschulStat *ka)
{
  gt_free(ka);
}

double gt_ka_get_lambda(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->lambda;
}

double gt_ka_get_logK(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->logK;
}

/* calculate probabilities of scores */
static ScoringFrequency *gt_karlin_altschul_stat_scoring_freuqnecy(
                                             const GtAlphabet *alphabet,
                                             const GtScoreHandler *scorehandler)
{
  unsigned int i, j, numofchars;
  GtWord score, range;

  ScoringFrequency *sf = gt_malloc(sizeof(*sf));
  gt_assert(sf);

  /* TODO: make generalizations of score/cost functionn,
   * min/max in scorehandler? */
  sf->low_score = -2;
  sf->high_score = 2;
  range = 2 - (-2) + 1; /*TODO: range = score_max-score_min+1*/
  sf->sprob = gt_malloc(range * sizeof (*sf->sprob));

  numofchars = gt_alphabet_num_of_chars(alphabet);

  for (i = 0; i < numofchars; i++)
  {
    for (j = 0; j < numofchars; j++)
      score = gt_scorehandler_get_replacement(scorehandler, i, j);

      if (score >= sf->low_score) /* TODO: error check */
        sf->sprob[score-sf->low_score] += nt_prob[i].p * nt_prob[j].p;
        /* TODO: make generalizations of alphabet probabilities,
           for now nt_prob */
  }

  return sf;
}

static double gt_karlin_altschul_stat_calculate_ungapped_lambda(
                                                     const ScoringFrequency *sf)
{
  double x0, x, lambda, tolerance, phi, dphi;
  GtWord i, low, high;

  /* solve phi(lambda) = -1 + sum_{i=l}^{u} sprob(i)*exp(i*lambda) = 0 */

  x0 = 0.5; /* x0 in (0,1) */
  tolerance = 1.e-5;
  GtWord kMaxIterations = 20;

  low = sf->low_score;
  high = sf->high_score;

  for (i = 0; i < kMaxIterations; i++)
  {
    phi = -1;
    for (i = low; i <= high; i++)
      phi += sf->sprob[i-low] * exp(i*x);

    dphi = 0;
    for (i = low; i <= high; i++)
      dphi += sf->sprob[i-low] * exp(i*x) * i;

    x = x0 - (phi)/(dphi);

    if (abs(x-x0) < tolerance)
      break;

    x0 = x;
  }

  lambda = x;

  /* TODO: better solution would be to apply newton's method to
   * polynomial (s. blast)*/

   return lambda;
}

static double gt_karlin_altschul_stat_calculate_H(const ScoringFrequency *sf,
                                                  double lambda)
{
  double H, sum, etonlambda;
  GtWord i, low, high;
  gt_assert(sf->sprob);

  low = sf->low_score;
  high = sf->high_score;

  etonlambda = exp(-lambda);
  sum = low * sf->sprob[low];
  for (i = low + 1; i <= high; i++)
    sum = i * sf->sprob[i-low] + etonlambda * sum;

  H = lambda * sum/pow(etonlambda,high);
  /*TODO: case underflow*/

   return H;
}

static double gt_karlin_altschul_stat_calculate_ungapped_K(const ScoringFrequency *sf,
                                                           double lambda)
{
  double H, K = 0;
  GtWord low, high, range, div;

  /* TODO: GT_Error object?
   * karlin-altschul theory works only if lambda > 0 && H > 0*/

  /* TODO: check score_average */

  H = gt_karlin_altschul_stat_calculate_H(sf, lambda);

  low = sf->low_score;
  high = sf->high_score;
  range = high - low;

  div = 1; /* TODO: greatest common divisor */

  low /= div;
  high /= div;
  lambda *= div;

  if (low == -1 && high == 1)
  {
    //TODO case 1
  }

  if (low == -1 || high == 1)
  {
    if (high != 1)
    {
      //Todo: case 3
    }

    //Todo: case 2
  }

  //TODO: otherwise case

   return K;
}

//TODO:new+fill oder trennen?
void gt_karlin_altschul_stat_calculate_params(GtKarlinAltschulStat *ka,
                                         GT_UNUSED bool ungapped_alignment,
                                         GT_UNUSED GtAlphabet *alphabet,
                                         GT_UNUSED GtScoreHandler *scorehandler)
{
  /* New ScoringFrequency */
  ScoringFrequency *sf =
                        gt_karlin_altschul_stat_scoring_freuqnecy(alphabet,
                                                                  scorehandler);

  /* karlin altschul parameters for ungapped alignments */
  ka->lambda = gt_karlin_altschul_stat_calculate_ungapped_lambda(sf);
  ka->K = gt_karlin_altschul_stat_calculate_ungapped_K(sf, ka->lambda);
  ka->logK = log(ka->K);
  ka->alpha_div_lambda = (1/ka->H);
  ka->beta = 0;

  /*TODO: gapped alignments*/
}
