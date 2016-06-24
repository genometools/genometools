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
#include "core/minmax.h"
#include "core/types_api.h"
#include "match/karlin_altschul_stat.h"
#include "extended/scorehandler.h"

/* TODO:reference, analog to blast */

struct GtKarlinAltschulStat
{
  double lambda,
         K,
         logK,
         H;
};

typedef struct{
  double *sprob,
         score_avg;
  GtWord low_score,
         high_score;
} ScoringFrequency;

typedef struct{
  char   ch;
  double p;
} LetterProb;

/* provisional solution */
static LetterProb nt_prob[] = {
  { 'A', 0.25 },
  { 'C', 0.25 },
  { 'G', 0.25 },
  { 'T', 0.25 }
};/* TODO: in generale normalize */

GtKarlinAltschulStat *gt_karlin_altschul_stat_new(void)
{
  GtKarlinAltschulStat *ka;
  ka = gt_malloc(sizeof (GtKarlinAltschulStat));
  ka->lambda = 0;
  ka->K = 0;
  ka->logK = 0;
  ka->H = 0;
  return ka;
}

void gt_karlin_altschul_stat_delete(GtKarlinAltschulStat *ka)
{
  gt_free(ka);
}

double gt_karlin_altschul_stat_get_lambda(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->lambda;
}

double gt_karlin_altschul_stat_get_logK(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->logK;
}

double gt_karlin_altschul_stat_get_K(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->K;
}

double gt_karlin_altschul_stat_get_alpha_div_lambda(const GtKarlinAltschulStat *ka,
                                                    GtWord matchscore,
                                                    GtWord mismatchscore)
{
  gt_assert(ka);
  if (matchscore == 1 && mismatchscore == 0)
    return 0;
  else
    return ka->lambda/ka->H;
  
  /* TODO: for ungapped */
}

double gt_karlin_altschul_stat_get_beta(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return 0; /* TODO: for ungapped */
}

/* calculate probabilities of scores */
static ScoringFrequency *gt_karlin_altschul_stat_scoring_frequency(
                                             const GtAlphabet *alphabet,
                                             const GtScoreHandler *scorehandler)
{
  unsigned int idx, jdx, numofchars;
  GtWord score, obs_min = 0, obs_max = 0, score_sum, range;
  double score_avg;
  
  gt_assert(alphabet && scorehandler);

  ScoringFrequency *sf = gt_malloc(sizeof(*sf));
  gt_assert(sf);

  numofchars = gt_alphabet_num_of_chars(alphabet);
  for (idx = 0; idx < numofchars; idx++)
  {
    for (jdx = 0; jdx < numofchars; jdx++)
    {
      score = gt_scorehandler_get_replacement(scorehandler, idx, jdx);
      obs_min = MIN(obs_min, score);
      obs_max = MAX(obs_max, score);
    }
  }

  /* for theoretically valid scoring systems */
  gt_assert(obs_min <= 0 && obs_max >= 0);
  sf->low_score = obs_min;
  sf->high_score = obs_max;
  
  range = obs_max - obs_min + 1;
  sf->sprob = gt_calloc(range, sizeof (*sf->sprob));

  for (idx = 0; idx < numofchars; idx++)
  {
    for (jdx = 0; jdx < numofchars; jdx++)
    {
      score = gt_scorehandler_get_replacement(scorehandler, idx, jdx);

      if (score >= sf->low_score)
        sf->sprob[score-sf->low_score] += nt_prob[idx].p * nt_prob[jdx].p;
        /* TODO: make generalizations of alphabet probabilities,
           for now nt_prob */
    }
  }

  score_sum = 0;
  for (score = obs_min; score <= obs_max; score++)
  {
    if (sf->sprob[score-obs_min] > 0)
      score_sum += score;
  }
  
  score_avg = 0;
  for (score = obs_min; score <= obs_max; score++)
  {
    sf->sprob[score-obs_min] /= score_sum;
    score_avg += score * sf->sprob[score-obs_min];
  }
  sf->score_avg = score_avg;

  return sf;
}

static double gt_karlin_altschul_stat_calculate_ungapped_lambda(
                                                     const ScoringFrequency *sf)
{
  double x0, x, lambda, tolerance, q, dq;
  GtWord low, high;
  GtUword kMaxIterations, idx, jdx;

  /* solve phi(lambda) = -1 + sum_{i=l}^{u} sprob(i)*exp(i*lambda) = 0 */

  x0 = 0.5; /* x0 in (0,1) */
  tolerance = 1.e-5;
  kMaxIterations = 20;

  low = sf->low_score;
  high = sf->high_score;

 /* write phi as phi(lambda) = exp(u*lambda) * q(exp(-lambda)) and solve the
  * polynomial q by apply newton's method
  *
  * q(x) = -x^u + sum_{k=0}^{u-l} sprob(u-k)* x^k
  */
  for (idx = 0; idx < kMaxIterations; idx++)
  {
    q = -pow(x0,high);
    for (jdx = 0; jdx <= high-low; jdx++)
      q += sf->sprob[high-low-jdx] * pow(x0,jdx);

    dq = -high*pow(x0,(high-1));
    for (jdx = 1; jdx <= high-low; jdx++)
      dq += sf->sprob[high-low-jdx] * jdx * pow(x0,jdx-1);

    x = x0 - (q/dq);

    if (fabs(x-x0) < tolerance)
      break;

    x0 = x;
  }

  lambda = -log(x);

  /* better solution would be to apply Horner's rule for evaluating a
     polynomial and its derivative (s. blast), but for the moment it works */

   return lambda;
}

static double gt_karlin_altschul_stat_calculate_H(const ScoringFrequency *sf,
                                                  double lambda)
{
  double H, sum, etonlambda;
  GtWord idx, low, high;
  gt_assert(sf->sprob);

  low = sf->low_score;
  high = sf->high_score;

  etonlambda = exp(-lambda);
  sum = low * sf->sprob[0];
  for (idx = low + 1; idx <= high; idx++)
    sum = idx * sf->sprob[idx-low] + etonlambda * sum;

  H = lambda * sum/pow(etonlambda,high);
  /*TODO: case underflow*/

   return H;
}

static double gt_karlin_altschul_stat_calculate_ungapped_K(const ScoringFrequency *sf,
                                                           double lambda, double H)
{
  GtWord low, high, div;
  double score_avg, score_avg_div, one_minus_expnlambda, K;

  gt_assert(lambda > 0 && H > 0);

  score_avg = sf->score_avg;
  gt_assert(score_avg < 0.0);

  div = 1; /* TODO: greatest common divisor */

  low = sf->low_score/div;
  high = sf->high_score/div;
  lambda *= div;
  
  /* range = high - low; */

  if (low == -1 && high == 1)
  {
      K = (sf->sprob[0] - sf->sprob[sf->high_score-sf->low_score]) *
          (sf->sprob[0] - sf->sprob[sf->high_score-sf->low_score]) / sf->sprob[0];
  }
  else if (low == -1 || high == 1)
  {
    one_minus_expnlambda = 1-exp(-lambda);
    if (high != 1)
    {
      score_avg_div = score_avg / div;
      K = lambda * one_minus_expnlambda / H *(score_avg_div * score_avg_div);
    }
    else
    {
      K = H/lambda * one_minus_expnlambda;
    }
  }
  else
  {
    //TODO: otherwise case
    
    /* K = lambda*exp(-2*sigma)/(H*(1-exp(-lambda)) */
    gt_assert(false); /* not implemented yet */
  }

  return K;
}

//TODO:new+fill oder trennen?
void gt_karlin_altschul_stat_calculate_params(GtKarlinAltschulStat *ka,
                                              GT_UNUSED bool ungapped_alignment,
                                              GtAlphabet *alphabet,
                                              GtScoreHandler *scorehandler)
{
  /* New ScoringFrequency */
  ScoringFrequency *sf =
                        gt_karlin_altschul_stat_scoring_frequency(alphabet,
                                                                  scorehandler);

  /* karlin altschul parameters for ungapped alignments */
  ka->lambda = gt_karlin_altschul_stat_calculate_ungapped_lambda(sf);
  ka->H = gt_karlin_altschul_stat_calculate_H(sf, ka->lambda);
  ka->K = gt_karlin_altschul_stat_calculate_ungapped_K(sf, ka->lambda, ka->H);
  ka->logK = log(ka->K);

  /*TODO: gapped alignments*/
}
