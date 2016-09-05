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
  -gt_evalue_calculate_searchspace
  -gt_evalue_calculate -> evalue
  -gt_karlin_altschul_stat_delete
 */

static GtUword gt_evalue_calculate_raw_score(const GtKarlinAltschulStat *ka,
                                             double bit_score)
{
  double raw_score, lambda, logK;
  gt_assert(ka);

  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);

  raw_score =  (bit_score * log(2) + logK)/lambda;
  return round(raw_score);
}

static GtUword gt_evalue_calculate_length_adjustment(GtUword query_length,
                                                     GtUword db_length,
                                                     GtUword num_of_db_seqs,
                                                     double alpha_div_lambda,
                                                     double beta,
                                                     double K,
                                                     double logK)
{
  unsigned int idx;
  const int kMaxIterations = 20;
  double len_min = 0, len_max, len_next, len, len_bar;
  GtUword length_adjustment;
  bool converged = false;
  double space, nNm;

  /* l_max is the largest nonnegative solution of
     K * (m - l) * (n - N * l) > MAX(m,n) */

  space = gt_safe_mult_ulong(db_length, query_length)
                                              - MAX(query_length, db_length)/K;
  if (space < 0)
    return 0; /* length_adjustnment = 0 */

  nNm = query_length * num_of_db_seqs + db_length;
  /* quadratic formula */
  len_max = 2 * space / (nNm + sqrt(nNm * nNm - 4 * num_of_db_seqs * space));

  len_next = 0;

  for (idx = 0; idx < kMaxIterations; idx++)
  {
    len = len_next;
    len_bar = beta + alpha_div_lambda *
              (logK + log((query_length-len)*(db_length-num_of_db_seqs*len)));
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

  length_adjustment = (GtUword) len_min; /* floor(fixed point ) */
  if (converged)
  {
    len = ceil(len_min);
    if (len <= len_max)
    {
      if (alpha_div_lambda * (log(K) + log((query_length-len) *
                             (db_length-num_of_db_seqs*len))) + beta >= len)
        length_adjustment = (GtUword) len;
    }
  }

  return length_adjustment;
}

/* deprecated */
GtUword gt_evalue_calculate_searchspace_dbencseq(const GtKarlinAltschulStat *ka,
                                                 const GtEncseq *dbencseq,
                                                 GtUword query_idx_length)
{
  GtUword total_length_of_db,
          num_of_db_seqs,
          actual_db_length,
          actual_query_length,
          effective_db_length,
          effective_query_length,
          length_adjustment;
  double alpha_div_lambda, beta, K, logK;

  gt_assert(ka);
  alpha_div_lambda = gt_karlin_altschul_stat_get_alphadlambda(ka);

  beta = gt_karlin_altschul_stat_get_beta(ka);
  K = gt_karlin_altschul_stat_get_K(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);

  /* db length */
  total_length_of_db = gt_encseq_total_length(dbencseq);
  num_of_db_seqs = gt_encseq_num_of_sequences(dbencseq);
  /* db length without seperators */
  actual_db_length = total_length_of_db - (num_of_db_seqs - 1);

  /* query length */
  actual_query_length = query_idx_length;

  length_adjustment = gt_evalue_calculate_length_adjustment(actual_query_length,
                                                            actual_db_length,
                                                            num_of_db_seqs,
                                                            alpha_div_lambda,
                                                            beta,
                                                            K, logK);

  effective_query_length = actual_query_length - length_adjustment;
  effective_db_length = actual_db_length -
                                      (num_of_db_seqs * length_adjustment);

  return gt_safe_mult_ulong(effective_query_length, effective_db_length);
}

GtUword gt_evalue_calculate_searchspace(const GtKarlinAltschulStat *ka,
                                        GtUword total_length_of_db,
                                        GtUword num_of_db_seqs,
                                        GtUword query_idx_length)
{
  GtUword actual_db_length,
          actual_query_length,
          effective_db_length,
          effective_query_length,
          length_adjustment;
  double alpha_div_lambda, beta, K, logK;

  gt_assert(ka);
  alpha_div_lambda = gt_karlin_altschul_stat_get_alphadlambda(ka);

  beta = gt_karlin_altschul_stat_get_beta(ka);
  K = gt_karlin_altschul_stat_get_K(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);

  /* db length without seperators */
  actual_db_length = total_length_of_db - (num_of_db_seqs - 1);

  /* query length */
  actual_query_length = query_idx_length;

  length_adjustment = gt_evalue_calculate_length_adjustment(actual_query_length,
                                                            actual_db_length,
                                                            num_of_db_seqs,
                                                            alpha_div_lambda,
                                                            beta,
                                                            K, logK);

  effective_query_length = actual_query_length - length_adjustment;
  effective_db_length = actual_db_length -
                                      (num_of_db_seqs * length_adjustment);

  return gt_safe_mult_ulong(effective_query_length, effective_db_length);
}

double gt_evalue_calculate_on_bitscore(const GtKarlinAltschulStat *ka,
                                       double bit_score,
                                       GtUword searchspace)
{
  double evalue, logK, lambda;
  GtUword raw_score;

  gt_assert(ka);
  raw_score = gt_evalue_calculate_raw_score(ka, bit_score);

  logK = gt_karlin_altschul_stat_get_logK(ka);
  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  evalue = searchspace * exp(-lambda*raw_score+logK);
  return evalue;
}

double gt_evalue_calculate(const GtKarlinAltschulStat *ka,
                           const GtScoreHandler *scorehandler,
                           GtUword ma,
                           GtUword mm,
                           GtUword id,
                           GtUword searchspace)
{
  double evalue, logK, lambda;
  GtWord matchscore, mismatchscore, gapscore;
  GtUword raw_score;

  gt_assert(ka && scorehandler);

  gt_assert(gt_scorehandler_get_gap_opening(scorehandler) == 0);
  /* only implemented for linear */

  matchscore = gt_scorehandler_get_matchscore(scorehandler);
  mismatchscore = gt_scorehandler_get_mismatchscore(scorehandler);
  gapscore = gt_scorehandler_get_gapscore(scorehandler);

  raw_score = ma*matchscore + mm*mismatchscore + id*gapscore;
  logK = gt_karlin_altschul_stat_get_logK(ka);
  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  evalue = searchspace*exp(-lambda*raw_score+logK);

  return evalue;
}

int gt_evalue_unit_test(GtError *err)
{
  GtKarlinAltschulStat *ka;
  GtScoreHandler *scorehandler;
  GtUword searchspace;
  double evalue_variance;

  int had_err = 0;
  gt_error_check(err);

  scorehandler = gt_scorehandler_new(1,-2,0,-2);
  ka = gt_karlin_altschul_stat_new(true, NULL, scorehandler, err);
  gt_error_check(err);

  /* checks searchspace calculation */
  gt_ensure(gt_evalue_calculate_searchspace(ka, 772376, 1952, 450)== 308243802);
  gt_ensure(gt_evalue_calculate_searchspace(ka, 772376, 1952, 300)== 199707252);
  gt_ensure(gt_evalue_calculate_searchspace(ka, 772376, 1952, 475)== 324731250);

  searchspace = gt_evalue_calculate_searchspace(ka, 772376, 1952, 300);

  /* checks evalue calculation */
  evalue_variance =
                  gt_evalue_calculate(ka, scorehandler, 300, 0, 0, searchspace)
                   /(6.148125*pow(10,-148));
  gt_ensure(evalue_variance > 0.99 && evalue_variance < 1.01);
  evalue_variance =
                  gt_evalue_calculate(ka, scorehandler, 213, 25, 1, searchspace)
                   /(4.220782*pow(10,-76));
  gt_ensure(evalue_variance > 0.99 && evalue_variance < 1.01);
  evalue_variance =
                  gt_evalue_calculate(ka, scorehandler, 206, 23, 1, searchspace)
                   /(1.499078*pow(10,-74));
  gt_ensure(evalue_variance > 0.99 && evalue_variance < 1.01);

  gt_scorehandler_delete(scorehandler);
  gt_karlin_altschul_stat_delete(ka);

  return had_err;
}
