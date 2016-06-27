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
#include "core/ma.h"
#include "core/minmax.h"
#include "core/safearith.h"
#include "core/types_api.h"
#include "match/evalue.h"
#include "match/karlin_altschul_stat.h"
#include "querymatch.h"

/* TODO:reference, analog to blast */

/*
 * information for invoking procedure:
 * new GtKarlinAltschulStat,
 * calculate GtKarlinAltschulStat params (1+2 together?)
 * -> calculate evalue
 * delete GtKarlinAltschulStat
 *
 */

static double gt_evalue_calculate_bit_score(const GtKarlinAltschulStat *ka,
                                            double raw_score)
{
  double bit_score, lambda, logK;
  gt_assert(ka);

  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);

  bit_score =  (raw_score * lambda - logK)/log(2);
  return bit_score;
}

/* TODO: reference to blast */
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

  
  space = gt_safe_mult_ulong(db_length, query_length) - MAX(query_length, db_length)/K;
  if (space < 0)
    return 0; /* length_adjustnment = 0 */

  nNm = query_length * num_of_db_seqs + db_length;
  len_max = 2 * space / (nNm + sqrt(nNm * nNm - 4 * num_of_db_seqs * space));  /* quadratic formula */

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

GtUword gt_evalue_calculate_searchspace(const GtEncseq *dbencseq,
                                        const GtEncseq *queryencseq,
                                        const GtKarlinAltschulStat *ka)
{
  GtUword total_length_of_db,
          total_length_of_query,
          num_of_db_seqs,
          num_of_query_seqs,
          actual_db_length,
          actual_query_length,
          effective_db_length,
          effective_query_length,
          length_adjustment;
  double alpha_div_lambda, beta, K, logK;

  gt_assert(ka);
  //TODO: berechnung, match, mismatch ergaenzen
  alpha_div_lambda = gt_karlin_altschul_stat_get_alpha_div_lambda(ka,1,0);
  /* 1,0 only useful for unit cost, TODO: generalize */

  beta = gt_karlin_altschul_stat_get_beta(ka);
  K = gt_karlin_altschul_stat_get_K(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);

  /* db length */
  total_length_of_db = gt_encseq_total_length(dbencseq);
  num_of_db_seqs = gt_encseq_num_of_sequences(dbencseq);
  /* db length without seperators */
  actual_db_length = total_length_of_db - (num_of_db_seqs - 1);

  /* query length */
  total_length_of_query = gt_encseq_total_length(queryencseq);
  num_of_query_seqs = gt_encseq_num_of_sequences(queryencseq);
  /* query length without seperators */
  actual_query_length = total_length_of_query - (num_of_query_seqs - 1);

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

double gt_evalue_calculate(const GtKarlinAltschulStat *ka,
                           const GtQuerymatch *querymatch,
                           GtUword searchspace)
{
  GtWord raw_score;
  double bit_score,
         pow2Nbitscore;
  
  gt_assert(ka);
  
  raw_score = gt_querymatch_score(querymatch);
  bit_score = gt_evalue_calculate_bit_score(ka, raw_score); 

  pow2Nbitscore = pow(2, -bit_score);
  gt_assert(pow2Nbitscore > DBL_MAX/searchspace);

  return searchspace * pow2Nbitscore;
}
