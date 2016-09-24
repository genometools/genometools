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
#include "core/unused_api.h"
#include "core/ensure.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "extended/scorehandler.h"
#include "match/karlin_altschul_stat.h"

#define K_ITER_MAX 100
#define K_SUMLIMIT_DEFAULT 0.0001

/*
  this library implements calculation of karlin-altschul parameter for E-value
  of Alignments analog to NCBI tool BLAST:

  Altschul S.F., Gish W., Miller W., Myers E.W. and Lipman D.J. (1990)
  Basic local alignment search tool. J. Mol. Biol. 215: 403-410.
 */

/* stores karlin altschul parameters */
struct GtKarlinAltschulStat
{
  double lambda,
         K,
         logK,
         H,
         alpha_div_lambda,
         beta;
  GtWord matchscore, mismatchscore, gapscore;
  GtUword total_length_db, num_of_db_seqs;
};

typedef struct{
  double *sprob,
         score_avg;
  GtWord low_align_score,
         high_align_score;
} ScoringFrequency;

/*
  precomputed values
  analog to BLAST

  1. Gap score,
  2. Lambda,
  3. K,
  4. H,
  5. Alpha,
  6. Beta,
 */

typedef enum
{
  gapextdidx,
  lambdaidx,
  Kidx,
  Hidx,
  alphaidx,
  betaidx
} GT_ValuesindeX;

/* the Blast Implementation uses matrices with more than one row, namely
   a row for each gap extension. If required, this should be adjusted
   accordingly. */

/* matchscore = 1 && mismatchscore = -4 */
static const double ga_vector_1_4[] = {
    -2,  1.26,  0.43, 0.90,  1.4, -1
};

/* matchscore = 2 && mismatchscore = -7 */
static const double ga_vector_2_7[] = {
    -4,  0.63, 0.43, 0.90,   0.7, -1
};

/* matchscore = 1 && mismatchscore = -3 */
static const double ga_vector_1_3[] = {
    -2,  1.25,  0.42, 0.83,  1.5, -2
};

/* matchscore = 2 && mismatchscore = -5 */
static const double ga_vector_2_5[] = {
    -4,  0.62, 0.39, 0.78,  0.8, -2
};

/* matchscore = 1 && mismatchscore = -2 */
static const double ga_vector_1_2[] = {
    -2, 1.19, 0.34, 0.66, 1.8, -3
};

/* matchscore = 2 && mismatchscore = -3 */
static const double ga_vector_2_3[] = {
    -4,  0.55, 0.21, 0.46,  1.2, -5
};

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

double gt_karlin_altschul_stat_get_alphadlambda(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->alpha_div_lambda;
}

double gt_karlin_altschul_stat_get_beta(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->beta;
}

/* calculate probabilities of scores */
static ScoringFrequency *gt_karlin_altschul_stat_scoring_frequency(
                                             unsigned int numofchars,
                                             const GtScoreHandler *scorehandler)
{
  unsigned int idx, jdx;
  GtWord score, obs_min = 0, obs_max = 0, range;
  double score_avg, score_sum;
  ScoringFrequency *sf;
  static double nt_prob[] = {
    0.25,
    0.25,
    0.25,
    0.25
  };

  gt_assert(scorehandler != NULL);
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
  sf = gt_malloc(sizeof (*sf));
  gt_assert(sf);
  sf->low_align_score = obs_min;
  sf->high_align_score = obs_max;

  range = obs_max - obs_min + 1;
  sf->sprob = gt_calloc(range, sizeof (*sf->sprob));

  for (idx = 0; idx < numofchars; idx++)
  {
    for (jdx = 0; jdx < numofchars; jdx++)
    {
      score = gt_scorehandler_get_replacement(scorehandler, idx, jdx);
      gt_assert(score >= obs_min);
      if (score >= obs_min)
      {
        sf->sprob[score - obs_min] += nt_prob[idx] * nt_prob[jdx];
      }
    }
  }

  score_sum = 0.0;
  for (score = obs_min; score <= obs_max; score++)
  {
    if (sf->sprob[score - obs_min] > 0)
    {
      score_sum += sf->sprob[score - obs_min];
    }
  }

  score_avg = 0.0;
  for (score = obs_min; score <= obs_max; score++)
  {
    sf->sprob[score - obs_min] /= score_sum;
    score_avg += score * sf->sprob[score - obs_min];
  }
  sf->score_avg = score_avg;
  return sf;
}

static double gt_karlin_altschul_stat_calculate_ungapped_lambda(
                                                     const ScoringFrequency *sf)
{
  double x0, x, lambda, tolerance, q, dq;
  GtWord low, high;
  GtUword k_max_iterations, idx, jdx;

  /* solve phi(lambda) = -1 + sum_{i=l}^{u} sprob(i)*exp(i*lambda) = 0 */

  x0 = 0.5; /* x0 in (0,1) */
  tolerance = 1.e-5;
  k_max_iterations = 20;

  low = sf->low_align_score;
  high = sf->high_align_score;

  /* write phi as phi(lambda) = exp(u*lambda) * q(exp(-lambda)) and solve the
    polynomial q by apply newton's method
    q(x) = -x^u + sum_{k=0}^{u-l} sprob(u-k)* x^k */
  for (idx = 0; idx < k_max_iterations; idx++)
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
     polynomial and its derivative (s. BLAST), but for the moment it works */

  return lambda;
}

static double gt_karlin_altschul_stat_calculate_H(const ScoringFrequency *sf,
                                                  double lambda)
{
  double H, sum, etonlambda;
  GtWord idx, low, high, scale;
  gt_assert(sf->sprob);

  low = sf->low_align_score;
  high = sf->high_align_score;

  etonlambda = exp(-lambda);
  sum = low * sf->sprob[0];
  for (idx = low + 1; idx <= high; idx++)
  {
    sum = idx * sf->sprob[idx-low] + etonlambda * sum;
  }

  scale = pow(etonlambda,high);
  if (scale > 0.0)
    H = lambda * sum/scale;
  else /* case underflow */
    H = lambda * exp(lambda * high + log(sum));

  return H;
}

static GtWord gt_karlin_altschul_stat_gcd(const ScoringFrequency *sf)
{
  GtWord idx, range, div;

  range = sf->high_align_score - sf->low_align_score+1;
  div = -sf->low_align_score;
  for (idx = 1; idx < range && div > 1; idx++)
  {
    if (sf->sprob[idx] != 0.0)
    {
      GtWord val = labs(idx+sf->low_align_score);
      if (val > div)
      {
        GtWord tmp = div;
        div = val;
        val = tmp;
      }
      while (val != 0)
      {
        GtWord tmp = div % val;
        div = val;
        val = tmp;
      }
    }
  }
  return div;
}

static double gt_karlin_altschul_stat_calculate_ungapped_K(
                                                     const ScoringFrequency *sf,
                                                     double lambda,
                                                     double H)
{
  GtWord  low,
          high,
          div,
          low_align_score,
          high_align_score,
          count, idx, jdx, firstidx, lastidx, secondidx,
          first, last;
  GtUword range,
          iterlimit,
          size,
          sigma = 0;
  double  score_avg,
          score_avg_div,
          one_minus_expnlambda,
          *alignnment_score_probs,
          expnlambda,
          K,
          sumlimit,
          inner_sum;

  gt_assert(lambda > 0 && H > 0);

  score_avg = sf->score_avg;
  gt_assert(score_avg < 0.0);

  /* greatest common divisor */
  div = gt_karlin_altschul_stat_gcd(sf);

  low = sf->low_align_score/div;
  high = sf->high_align_score/div;
  lambda *= div;

  range = high - low;
  expnlambda = exp(-lambda);

  if (low == -1 && high == 1)
  {
    K = (sf->sprob[0] - sf->sprob[sf->high_align_score-sf->low_align_score]) *
        (sf->sprob[0] - sf->sprob[sf->high_align_score-sf->low_align_score])/
        sf->sprob[0];
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
    /* K = lambda*exp(-2*sigma)/(H*(1-exp(-lambda)) */

    sumlimit = K_SUMLIMIT_DEFAULT;
    iterlimit = K_ITER_MAX;

    size = iterlimit * range + 1;
    alignnment_score_probs = gt_calloc(size, sizeof (*alignnment_score_probs));
    gt_assert(alignnment_score_probs);

    low_align_score = 0;
    high_align_score = 0;
    inner_sum = 1.0;
    alignnment_score_probs[0] = 1.0;

    for (count = 0; count < iterlimit && inner_sum > sumlimit; count++)
    {
      if (count > 0)
      {
        inner_sum /= count;
        sigma += inner_sum;
      }

      first = last = range;
      low_align_score += low;
      high_align_score += high;
      for (idx = high_align_score-low_align_score; idx >= 0; idx--)
      {
        firstidx = idx-first;
        lastidx = idx-last;
        secondidx = sf->sprob[low] + first;
        for (inner_sum = 0.; firstidx >= lastidx; )
        {
          inner_sum += alignnment_score_probs[firstidx] *
                      alignnment_score_probs[secondidx];
          firstidx--;
          secondidx++;
        }

        if (first > 0)
          --first;
        if (idx <= range)
          --last;

        alignnment_score_probs[idx] = inner_sum;
        if (idx == 0)
          break;
      }

      inner_sum = alignnment_score_probs[++idx];
      for (jdx = low_align_score + 1; jdx < 0; jdx++)
      {
        inner_sum = alignnment_score_probs[++idx] + inner_sum * expnlambda;
      }
      inner_sum *= expnlambda;
      for (/*Nothing*/; jdx <= high_align_score; ++jdx)
        inner_sum += alignnment_score_probs[++jdx];
    }
    gt_free(alignnment_score_probs);
    /* no terms of geometric progression, check to add these terms for
       correction in future */

    K = -exp(-2.0*sigma)/(H/lambda*expnlambda);
  }

  return K;
}

static void get_values_from_vector(GtKarlinAltschulStat *ka,
                                   const double *vector)
{
  ka->lambda = vector[lambdaidx];
  ka->K = vector[Kidx];
  ka->logK = log(ka->K);
  ka->H = vector[Hidx];
  gt_assert(ka->lambda != 0.0);
  ka->alpha_div_lambda = vector[alphaidx]/ka->lambda;
  ka->beta = vector[betaidx];
}

static void gt_karlin_altschul_stat_get_gapped_params(GtKarlinAltschulStat *ka)
{
  const double *ga_vector = NULL;

  gt_assert(ka != NULL);
  if (ka->matchscore == 1 && ka->mismatchscore == -4)
  {
    ga_vector = ga_vector_1_4;
  }
  else if (ka->matchscore == 2 && ka->mismatchscore == -7)
  {
    ga_vector = ga_vector_2_7;
  }
  else if (ka->matchscore == 1 && ka->mismatchscore == -3)
  {
    ga_vector = ga_vector_1_3;
  }
  else if (ka->matchscore == 2 && ka->mismatchscore == -5)
  {
    ga_vector = ga_vector_2_5;
  }
  else if (ka->matchscore == 1 && ka->mismatchscore == -2)
  {
    ga_vector = ga_vector_1_2;
  }
  else if (ka->matchscore == 2 && ka->mismatchscore == -3)
  {
    ga_vector = ga_vector_2_3;
  }
  else
  {
    fprintf(stderr,"no precomputed values for combination matchscore "
                   GT_WD " and mismatchscore " GT_WD " in evalue calculation "
                   "of gapped alignments",ka->matchscore, ka->mismatchscore);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_assert(ga_vector[gapextdidx] == ka->gapscore);
  get_values_from_vector(ka, ga_vector);
}

GtKarlinAltschulStat *gt_karlin_altschul_stat_new(unsigned int numofchars,
                                                  const GtScoreHandler
                                                        *scorehandler)
{
  GtKarlinAltschulStat *ka = gt_malloc(sizeof *ka);

  gt_assert(ka);
  ka->lambda = 0;
  ka->K = 0;
  ka->logK = 0;
  ka->H = 0;
  ka->total_length_db = GT_UWORD_MAX;
  ka->num_of_db_seqs = GT_UWORD_MAX;
  gt_assert(gt_scorehandler_get_gap_opening(scorehandler) == 0);
  ka->matchscore = gt_scorehandler_get_matchscore(scorehandler);
  ka->mismatchscore = gt_scorehandler_get_mismatchscore(scorehandler);
  ka->gapscore = gt_scorehandler_get_gapscore(scorehandler);
  /* only implemented for linear scores */
  if (numofchars == 0)
  { /* gapped alignment */
    gt_karlin_altschul_stat_get_gapped_params(ka);
  } else
  {
    /* New ScoringFrequency */
    ScoringFrequency *sf
      = gt_karlin_altschul_stat_scoring_frequency(numofchars,scorehandler);
    gt_assert(sf != NULL);

    /* karlin altschul parameters for ungapped alignments */
    ka->lambda = gt_karlin_altschul_stat_calculate_ungapped_lambda(sf);
    ka->H = gt_karlin_altschul_stat_calculate_H(sf, ka->lambda);
    ka->K = gt_karlin_altschul_stat_calculate_ungapped_K(sf, ka->lambda, ka->H);
    ka->logK = log(ka->K);
    gt_assert(ka->H != 0.0);
    ka->alpha_div_lambda = 1/ka->H;
    ka->beta = 0;
    gt_free(sf->sprob);
    gt_free(sf);
  }
  return ka;
}

GtKarlinAltschulStat *gt_karlin_altschul_stat_new_gapped(void)
{
  const unsigned int gapped_alignment_flag = 0;
  GtScoreHandler *scorehandler = gt_scorehandler_new(1,-2,0,-2);
  GtKarlinAltschulStat *karlin_alt_schul_stat
    = gt_karlin_altschul_stat_new(gapped_alignment_flag,scorehandler);
  gt_scorehandler_delete(scorehandler);
  return karlin_alt_schul_stat;
}

void gt_karlin_altschul_stat_add_keyvalues(
                             GtKarlinAltschulStat *karlin_altschul_stat,
                             GtUword total_length_db,
                             GtUword num_of_db_seqs)
{
  if (karlin_altschul_stat != NULL)
  {
    gt_assert(karlin_altschul_stat->total_length_db == GT_UWORD_MAX);
    karlin_altschul_stat->total_length_db = total_length_db;
    karlin_altschul_stat->num_of_db_seqs = num_of_db_seqs;
  }
}

GtWord gt_karlin_altschul_stat_mismatchscore(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->mismatchscore;
}

GtWord gt_karlin_altschul_stat_matchscore(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->matchscore;
}

GtWord gt_karlin_altschul_stat_gapscore(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->gapscore;
}

GtWord gt_karlin_altschul_get_total_length_db(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL && ka->total_length_db != GT_UWORD_MAX);
  return ka->total_length_db;
}

GtWord gt_karlin_altschul_get_num_of_db_seqs(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL && ka->num_of_db_seqs != GT_UWORD_MAX);
  return ka->num_of_db_seqs;
}

int gt_karlin_altschul_stat_unit_test(GtError *err)
{
  GtKarlinAltschulStat *ka;
  GtScoreHandler *scorehandler;
  double q;
  unsigned int numofchars = 0;

  int had_err = 0;
  gt_error_check(err);

  scorehandler = gt_scorehandler_new(1,-2,0,-2);

  /* check function for gapped alignments */
  ka = gt_karlin_altschul_stat_new(0, scorehandler);
  gt_ensure(ka->lambda == 1.19);
  gt_ensure(ka->H == 0.66);
  gt_ensure(ka->K == 0.34);
  gt_karlin_altschul_stat_delete(ka);

  /* check function for ungapped alignments */
  numofchars = 4;
  ka = gt_karlin_altschul_stat_new(numofchars, scorehandler);
  q = ka->lambda/1.33; /* lambda = 1.33 */
  gt_ensure(0.99 < q && q < 1.01);
  q = ka->H/1.12; /* H = 1.12 */
  gt_ensure(0.99 < q && q < 1.01);
  q = ka->K/0.621; /* K = 0.621 */
  gt_ensure(0.99 < q && q < 1.01);
  gt_karlin_altschul_stat_delete(ka);

  gt_scorehandler_delete(scorehandler);
  return had_err;
}
