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
#include "core/ensure.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/safearith.h"
#include "core/types_api.h"
#include "match/evalue.h"
#include "match/karlin_altschul_stat.h"

/*
  this library implements calculation of E-value of Alignments analog to NCBI
  tool BLAST:

    Altschul S.F., Gish W., Miller W., Myers E.W. and Lipman D.J. (1990)
    Basic local alignment search tool. J. Mol. Biol. 215: 403-410.
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
     K * (m - l) * (n - N * l) > MAX(m,n) */

  space = gt_safe_mult_ulong(actual_db_length, query_length)
          - MAX(query_length, actual_db_length)/K;
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
              (logK + log((query_length-len)*
                          (actual_db_length-num_of_db_seqs*len)));
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
      if (alpha_div_lambda * (log(K) +
          log((query_length-len) * (actual_db_length-num_of_db_seqs*len)))
          + beta >= len)
      {
        length_adjustment = (GtUword) len;
      }
    }
  }

  return length_adjustment;
}

GtUword gt_evalue_searchspace(const GtKarlinAltschulStat *ka,
                              GtUword query_idx_length)
{
  GtUword actual_db_length,
          num_of_db_seqs,
          actual_query_length,
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
  gt_assert(num_of_db_seqs > 0);
  actual_db_length = gt_karlin_altschul_get_total_length_db(ka) -
                     (num_of_db_seqs - 1);

  /* query length */
  actual_query_length = query_idx_length;

  length_adjustment = gt_evalue_length_adjustment(actual_query_length,
                                                  actual_db_length,
                                                  num_of_db_seqs,
                                                  alpha_div_lambda,
                                                  beta,
                                                  K,
                                                  logK);

  effective_query_length = actual_query_length - length_adjustment;
  effective_db_length
    = actual_db_length - (num_of_db_seqs * length_adjustment);

  return gt_safe_mult_ulong(effective_query_length, effective_db_length);
}

GtWord gt_evalue_raw_score(const GtKarlinAltschulStat *ka,
                           GtUword matches,
                           GtUword mismatches,
                           GtUword indels)
{
  GtWord matchscore, mismatchscore, gapscore;

  gt_assert(ka);
  matchscore = gt_karlin_altschul_stat_matchscore(ka);
  mismatchscore = gt_karlin_altschul_stat_mismatchscore(ka);
  gapscore = gt_karlin_altschul_stat_gapscore(ka);
  gt_assert(matchscore > 0 && mismatchscore < 0 && gapscore < 0);
  return (GtWord) matches * matchscore + (GtWord) mismatches * mismatchscore +
         (GtWord) indels * gapscore;
}

double gt_evalue_raw_score2bit_score(const GtKarlinAltschulStat *ka,
                                     GtWord raw_score)
{
  double lambda, logK;

  gt_assert(ka);
  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);
  return (lambda * (double) raw_score - logK)/log(2);
}

GtWord gt_evalue_bit_score2raw_score(const GtKarlinAltschulStat *ka,
                                     double bit_score)
{
  double raw_score, lambda, logK;

  gt_assert(ka);
  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);

  raw_score = (bit_score * log(2) + logK)/lambda;
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

  gt_assert(ka);
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

  gt_assert(ka);
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
  ka = gt_karlin_altschul_stat_new(numchars,scorehandler);
  gt_karlin_altschul_stat_add_keyvalues(ka,772376, 1952);

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
