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

#include <float.h>
#include <math.h>
#include "core/unused_api.h"
#include "core/ensure_api.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/minmax_api.h"
#include "core/types_api.h"
#include "core/encseq.h"
#include "core/radix_sort.h"
#include "extended/scorehandler.h"
#include "karlin_altschul_stat.h"

/*
  the first part of this file implements calculation of karlin-altschul
  parameter for E-value of Alignments in analogy to NCBI tool BLAST:

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
         beta,
         log2;
  GtWord matchscore, mismatchscore, gapscore;
  GtUword actual_length_db, num_of_db_seqs, num_of_query_seqs,
          different_lengths;
  GtUwordPair *searchspace_store;
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

/* matchscore = 1 && mismatchscore = -1 */
static const double ga_vector_1_1[] = {
   -2, 0.80, 0.064, 0.17, 4.8, -16
};

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
  if (ka != NULL)
  {
    gt_free(ka->searchspace_store);
    gt_free(ka);
  }
}

static double gt_karlin_altschul_stat_get_lambda(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->lambda;
}

static double gt_karlin_altschul_stat_get_logK(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->logK;
}

static double gt_karlin_altschul_stat_get_K(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->K;
}

static double gt_karlin_altschul_stat_get_alphadlambda(
                         const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->alpha_div_lambda;
}

static double gt_karlin_altschul_stat_get_beta(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
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
      obs_min = GT_MIN(obs_min, score);
      obs_max = GT_MAX(obs_max, score);
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
  GtUword k_max_iterations = 20, idx, jdx;

  /* solve phi(lambda) = -1 + sum_{i=l}^{u} sprob(i)*exp(i*lambda) = 0 */

  x0 = 0.5; /* x0 in (0,1) */
  tolerance = 1.e-5;

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

    dq = -high*pow(x0,high-1);
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
  GtWord idx, range, divisor;

  range = sf->high_align_score - sf->low_align_score+1;
  divisor = -sf->low_align_score;
  for (idx = 1; idx < range && divisor > 1; idx++)
  {
    if (sf->sprob[idx] != 0.0)
    {
      GtWord val = labs(idx+sf->low_align_score);
      if (val > divisor)
      {
        GtWord tmp = divisor;
        divisor = val;
        val = tmp;
      }
      while (val != 0)
      {
        GtWord tmp = divisor % val;
        divisor = val;
        val = tmp;
      }
    }
  }
  return divisor;
}

static double gt_karlin_altschul_stat_calculate_ungapped_K(
                                                     const ScoringFrequency *sf,
                                                     double lambda,
                                                     double H)
{
  GtWord  low,
          high,
          divisor;
  double  expnlambda,
          K;

  gt_assert(lambda > 0 && H > 0);
  divisor = gt_karlin_altschul_stat_gcd(sf); /* greatest common divisor */
  low = sf->low_align_score/divisor;
  high = sf->high_align_score/divisor;
  lambda *= divisor;
  expnlambda = exp(-lambda);
  if (low == -1 && high == 1)
  {
    K = (sf->sprob[0] - sf->sprob[sf->high_align_score - sf->low_align_score]) *
        (sf->sprob[0] - sf->sprob[sf->high_align_score - sf->low_align_score])/
        sf->sprob[0];
  }
  else if (low == -1 || high == 1)
  {
    const double score_avg = sf->score_avg,
                 one_minus_expnlambda = 1-expnlambda;

    gt_assert(score_avg < 0.0);
    if (high != 1)
    {
      const double score_avg_div = score_avg / divisor;
      K = lambda * one_minus_expnlambda / H *(score_avg_div * score_avg_div);
    }
    else
    {
      K = H/lambda * one_minus_expnlambda;
    }
  }
  else
  {
    const GtUword range = (GtUword) (high - low);
    const GtUword iterlimit = 100;
    const double sumlimit = 0.0001;
    const GtUword size = iterlimit * range + 1;
    GtWord count;
    GtUword sigma = 0;
    GtWord low_align_score = 0, high_align_score = 0;
    double inner_sum = 1.0,
           *alignnment_score_probs
             = gt_calloc(size, sizeof (*alignnment_score_probs));

    gt_assert(low <= high);
    /* K = lambda*exp(-2*sigma)/(H*(1-exp(-lambda)) */

    alignnment_score_probs[0] = 1.0;
    for (count = 0; count < iterlimit && inner_sum > sumlimit; count++)
    {
      GtWord idx, jdx, first, last;
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
        GtWord firstidx = idx-first;
        GtWord lastidx = idx-last;
        GtWord secondidx = sf->sprob[low] + first;
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
  gt_assert(ka->lambda != 0.0);
  ka->K = vector[Kidx];
  ka->logK = log(ka->K);
  ka->H = vector[Hidx];
  ka->alpha_div_lambda = vector[alphaidx]/ka->lambda;
  ka->beta = vector[betaidx];
}

static const double *gt_karlin_altschul_score_vector(GtWord matchscore,
                                                     GtWord mismatchscore)
{
  if (matchscore == 1 && mismatchscore == -1)
  {
    return ga_vector_1_1;
  }
  if (matchscore == 1 && mismatchscore == -4)
  {
    return ga_vector_1_4;
  }
  if (matchscore == 2 && mismatchscore == -7)
  {
    return ga_vector_2_7;
  }
  if (matchscore == 1 && mismatchscore == -3)
  {
    return ga_vector_1_3;
  }
  if (matchscore == 2 && mismatchscore == -5)
  {
    return ga_vector_2_5;
  }
  if (matchscore == 1 && mismatchscore == -2)
  {
    return ga_vector_1_2;
  }
  if (matchscore == 2 && mismatchscore == -3)
  {
    return ga_vector_2_3;
  }
  return NULL;
}

static void gt_karlin_altschul_stat_get_gapped_params(GtKarlinAltschulStat *ka)
{
  const double *ga_vector = NULL;

  gt_assert(ka != NULL);
  ga_vector = gt_karlin_altschul_score_vector(ka->matchscore,ka->mismatchscore);
  if (ga_vector == NULL)
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

  gt_assert(ka != NULL);
  ka->searchspace_store = NULL;
  ka->num_of_query_seqs = 0;
  ka->different_lengths = 0;
  ka->lambda = 0;
  ka->K = 0;
  ka->logK = 0;
  ka->H = 0;
  ka->log2 = log(2);
  ka->actual_length_db = GT_UWORD_MAX;
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

GtKarlinAltschulStat *gt_karlin_altschul_stat_new_gapped(
                             GtUword total_length_db,
                             GtUword num_of_db_seqs,
                             GT_UNUSED const GtEncseq *query_encseq)
{
  const unsigned int gapped_alignment_flag = 0;
  GtScoreHandler *scorehandler = gt_scorehandler_new(1,-2,0,-2);
  GtKarlinAltschulStat *ka
    = gt_karlin_altschul_stat_new(gapped_alignment_flag,scorehandler);

  gt_scorehandler_delete(scorehandler);
  gt_assert(num_of_db_seqs > 0 &&
            ka->actual_length_db == GT_UWORD_MAX);
  ka->actual_length_db = total_length_db - (num_of_db_seqs - 1);
  ka->num_of_db_seqs = num_of_db_seqs;
  if (query_encseq != NULL)
  {
    GtUword *seq_length_tab = gt_all_sequence_lengths_get(query_encseq);
    GtUwordPair *store_under_construction;

    ka->num_of_query_seqs = gt_encseq_num_of_sequences(query_encseq);
    if (seq_length_tab == NULL) /* all are of the same length */
    {
      store_under_construction = gt_malloc(sizeof *ka->searchspace_store * 1);
      store_under_construction[0].a = gt_encseq_seqlength(query_encseq,0);
      store_under_construction[0].b
        = gt_evalue_searchspace(ka,store_under_construction[0].a);
      ka->different_lengths = 1;
    } else
    {
      GtUword ridx, widx;

      gt_radixsort_inplace_ulong(seq_length_tab,ka->num_of_query_seqs);
      gt_assert(ka->num_of_query_seqs > 0);
      for (ridx = 1, widx = 0; ridx < ka->num_of_query_seqs; ridx++)
      {
        if (seq_length_tab[ridx] != seq_length_tab[widx])
        {
          seq_length_tab[++widx] = seq_length_tab[ridx];
        }
      }
      ka->different_lengths = widx + 1;
      store_under_construction = gt_malloc(sizeof *ka->searchspace_store *
                                           ka->different_lengths);
      for (widx = 0; widx < ka->different_lengths; widx++)
      {
        store_under_construction[widx].a = seq_length_tab[widx];
        store_under_construction[widx].b
          = gt_evalue_searchspace(ka,store_under_construction[widx].a);
      }
      gt_free(seq_length_tab);
    }
    ka->searchspace_store = store_under_construction;
  }
  return ka;
}

static GtWord gt_karlin_altschul_stat_mismatchscore(
                                       const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->mismatchscore;
}

static GtWord gt_karlin_altschul_stat_matchscore(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->matchscore;
}

static GtWord gt_karlin_altschul_stat_gapscore(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL);
  return ka->gapscore;
}

static GtWord gt_karlin_altschul_get_actual_length_db(
                       const GtKarlinAltschulStat *ka)
{
  gt_assert(ka != NULL && ka->actual_length_db != GT_UWORD_MAX);
  return ka->actual_length_db;
}

static GtWord gt_karlin_altschul_get_num_of_db_seqs(
               const GtKarlinAltschulStat *ka)
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
  ka = gt_karlin_altschul_stat_new(0, scorehandler); /* unit test */
  gt_ensure(ka->lambda == 1.19);
  gt_ensure(ka->H == 0.66);
  gt_ensure(ka->K == 0.34);
  gt_karlin_altschul_stat_delete(ka);

  /* check function for ungapped alignments */
  numofchars = 4;
  ka = gt_karlin_altschul_stat_new(numofchars, scorehandler); /* unit test */
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

/*
  the rest of this file implements function for calculation of E-value of
  Alignments.
 */

/*
  information for invoking procedure:
  -gt_karlin_altschul_stat_new
  -gt_evalue_searchspace
  -gt_evalue_from_eop_count -> evalue
  -gt_karlin_altschul_stat_delete
 */

static GtUword gt_evalue_length_adjustment(GtUword query_length,
                                           GtUword actual_db_length,
                                           GtUword num_of_db_seqs,
                                           double alpha_div_lambda,
                                           double beta,
                                           double K,
                                           double logK)
{
  unsigned int idx;
  const unsigned int kMaxIterations = 20;
  double len_min = 0, len_max, len_next, len, len_bar, space, nNm;
  GtUword length_adjustment;
  bool converged = false;

  /* l_max is the largest nonnegative solution of
     K * (m - l) * (n - N * l) > GT_MAX(m,n) */

  space = actual_db_length * query_length
          - GT_MAX(query_length, actual_db_length)/K;
  if (space < 0)
    return 0; /* length_adjustnment = 0 */

  nNm = query_length * num_of_db_seqs + actual_db_length;
  /* quadratic formula */
  len_max = 2 * space / (nNm + sqrt(nNm * nNm - 4 * num_of_db_seqs * space));

  len_next = 0;

  for (idx = 0; idx < kMaxIterations; idx++)
  {
    len = len_next;
    len_bar = beta + alpha_div_lambda *
              (logK + log((double) (query_length - len) *
                          (double) (actual_db_length - num_of_db_seqs*len)));
    if (len_bar >= len)
    {
      len_min = len;
      if (len_bar -len_min <= 1.0)
      {
        converged = true;
        break;
      }
      if (len_min == len_max)
        break;
    }
    else
    {
      len_max = len;
    }
    if (len_min <= len_bar && len_bar <= len_max)
      len_next = len_bar;
    else if (idx == 0)
      len_next = len_max;
    else
      len_next = (len_min+len_max)/2;
  }

  length_adjustment = (GtUword) len_min; /* floor(fixed point) */
  if (converged)
  {
    len = ceil(len_min);
    if (len <= len_max)
    {
      if (alpha_div_lambda * (logK +
                              log((query_length-len) *
                                  (actual_db_length-num_of_db_seqs*len))) + beta
          >= len)
      {
        length_adjustment = (GtUword) len;
      }
    }
  }

  return length_adjustment;
}

static GtUword gt_searchspace_store_find(const GtUwordPair *arr,
                                         GtUword len,
                                         GtUword key)
{
  const GtUwordPair *left = arr, *right = arr + len - 1;

  while (left <= right)
  {
    const GtUwordPair *mid = left + GT_DIV2(right - left);

    if (key < mid->a)
    {
      right = mid - 1;
    } else
    {
      if (key > mid->a)
      {
        left = mid + 1;
      } else
      {
        return mid->b;
      }
    }
  }
  fprintf(stderr,"%s: cannot find entry for " GT_WU "\n",__func__,key);
  gt_assert(false);
  return GT_UWORD_MAX;
}

GtUword gt_evalue_searchspace(const GtKarlinAltschulStat *ka,
                              GtUword query_idx_length)
{
  if (ka->searchspace_store != NULL)
  {
    return gt_searchspace_store_find(ka->searchspace_store,
                                     ka->different_lengths,
                                     query_idx_length);
  } else
  {
    GtUword actual_db_length,
            num_of_db_seqs,
            effective_db_length,
            effective_query_length,
            length_adjustment;
    double alpha_div_lambda, beta, K, logK;

    gt_assert(ka != NULL);
    alpha_div_lambda = gt_karlin_altschul_stat_get_alphadlambda(ka);

    beta = gt_karlin_altschul_stat_get_beta(ka);
    K = gt_karlin_altschul_stat_get_K(ka);
    logK = gt_karlin_altschul_stat_get_logK(ka);

    num_of_db_seqs = gt_karlin_altschul_get_num_of_db_seqs(ka);
    actual_db_length = gt_karlin_altschul_get_actual_length_db(ka);
    length_adjustment = gt_evalue_length_adjustment(query_idx_length,
                                                    actual_db_length,
                                                    num_of_db_seqs,
                                                    alpha_div_lambda,
                                                    beta,
                                                    K,
                                                    logK);

    effective_query_length = query_idx_length - length_adjustment;
    effective_db_length
      = actual_db_length - (num_of_db_seqs * length_adjustment);

    return effective_query_length * effective_db_length;
  }
}

GtWord gt_evalue_raw_score(const GtKarlinAltschulStat *ka,
                           GtUword matches,
                           GtUword mismatches,
                           GtUword indels)
{
  GtWord matchscore, mismatchscore, gapscore;

  gt_assert(ka != NULL);
  matchscore = gt_karlin_altschul_stat_matchscore(ka);
  mismatchscore = gt_karlin_altschul_stat_mismatchscore(ka);
  gapscore = gt_karlin_altschul_stat_gapscore(ka);
  gt_assert(matchscore > 0 && mismatchscore < 0 && gapscore < 0);
  return (GtWord) matches * matchscore +
         (GtWord) mismatches * mismatchscore +
         (GtWord) indels * gapscore;
}

double gt_evalue_raw_score2bit_score(const GtKarlinAltschulStat *ka,
                                     GtWord raw_score)
{
  double lambda, logK;

  gt_assert(ka != NULL);
  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);
  return (lambda * (double) raw_score - logK)/ka->log2;
}

GtWord gt_evalue_bit_score2raw_score(const GtKarlinAltschulStat *ka,
                                     double bit_score)
{
  double raw_score, lambda, logK;

  gt_assert(ka != NULL);
  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);

  raw_score = (bit_score * ka->log2 + logK)/lambda;
  return (GtWord) round(raw_score);
}

double gt_evalue_from_raw_score(const GtKarlinAltschulStat *ka,
                                GtWord raw_score,
                                GtUword searchspace)
{
  double logK = gt_karlin_altschul_stat_get_logK(ka);
  double lambda = gt_karlin_altschul_stat_get_lambda(ka);
  return searchspace * exp(-lambda * raw_score + logK);
}

double gt_evalue_from_bitscore(const GtKarlinAltschulStat *ka,
                               double bit_score,
                               GtUword searchspace)
{
  GtWord raw_score;

  gt_assert(ka != NULL);
  raw_score = gt_evalue_bit_score2raw_score(ka, bit_score);
  return gt_evalue_from_raw_score(ka,raw_score,searchspace);
}

double gt_evalue_from_eop_count(const GtKarlinAltschulStat *ka,
                                GtUword matches,
                                GtUword mismatches,
                                GtUword indels,
                                GtUword searchspace)
{
  GtWord raw_score;

  gt_assert(ka != NULL);
  raw_score = gt_evalue_raw_score(ka, matches, mismatches, indels);
  return gt_evalue_from_raw_score(ka,raw_score,searchspace);
}

int gt_evalue_unit_test(GT_UNUSED GtError *err)
{
  GtKarlinAltschulStat *ka;
  GtScoreHandler *scorehandler;
  GtUword searchspace;
  double evalue_variance;
  const unsigned int numchars = 0; /* means gapped case */
  int had_err = 0;

  scorehandler = gt_scorehandler_new(1,-2,0,-2);
  ka = gt_karlin_altschul_stat_new(numchars,scorehandler); /* unit test */
  ka->actual_length_db = 772376 - (1952 - 1);
  ka->num_of_db_seqs = 1952;

  /* checks searchspace calculation */
  gt_ensure(gt_evalue_searchspace(ka, 450)== 308243802);
  gt_ensure(gt_evalue_searchspace(ka, 300)== 199707252);
  gt_ensure(gt_evalue_searchspace(ka, 475)== 324731250);

  searchspace = gt_evalue_searchspace(ka, 300);

  /* checks evalue calculation */
  evalue_variance
    = gt_evalue_from_eop_count(ka, 300, 0, 0, searchspace)
      /(6.148125 * pow(10,-148));
  gt_ensure(evalue_variance > 0.99 && evalue_variance < 1.01);
  evalue_variance
    = gt_evalue_from_eop_count(ka, 213, 25, 1, searchspace)
      /(4.220782 * pow(10,-76));
  gt_ensure(evalue_variance > 0.99 && evalue_variance < 1.01);
  evalue_variance
    = gt_evalue_from_eop_count(ka, 206, 23, 1, searchspace)
      /(1.499078 * pow(10,-74));
  gt_ensure(evalue_variance > 0.99 && evalue_variance < 1.01);

  gt_scorehandler_delete(scorehandler);
  gt_karlin_altschul_stat_delete(ka);
  return had_err;
}
