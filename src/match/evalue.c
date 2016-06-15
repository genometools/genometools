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
 * calculate bitscore
 * calculate searchspace
 * -> calculate evalue
 * delete GtKarlinAltschulStat
 *
 */

double gt_evalue_calculate_bit_score(double raw_score,
                                     const GtKarlinAltschulStat *ka)
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
                                                     double K)
{
  unsigned int idx;
  const int kMaxIterations = 20;
  double l_min = 0, l_max, l_next, l, l_bar;
  GtUword length_adjustment;
  bool converged = false;

  /* l_max is the largest nonnegative solution of
     K * (m - l) * (n - N * l) > MAX(m,n) */

  
  double c = gt_safe_mult_ulong(db_length, query_length) - MAX(query_length, db_length)/K;
  if (c < 0)
    return 0; /* length_adjustnment = 0 */

  double a = num_of_db_seqs;
  double mb = query_length * num_of_db_seqs + db_length;
  l_max = 2 * c / (mb + sqrt(mb * mb - 4 * a * c));  /* quadratic formula */

  l_next = 0;
  for (idx = 0; idx < kMaxIterations; idx++)
  {
    l = l_next;
    l_bar = beta + alpha_div_lambda *
            (log(K) + log((query_length-l)*(db_length-num_of_db_seqs*l)));
    if (l_bar >= l)
    {
      l_min = l;
      if (l_bar -l_min <= 1.0)
      {
        converged = true;
        break;
      }
      if (l_min == l_max)
        break;
    }
    else
    {
      l_max = l;
    }
    if (l_min <= l_bar && l_bar <= l_max)
      l_next = l_bar;
    else if (idx == 0)
      l_next = l_max;
    else
      l_next = (l_min+l_max)/2;
  }

  length_adjustment = (GtUword) l_min; /* floor(fixed point ) */
  if (converged)
  {
    l = ceil(l_min);
    if (l <= l_max)
    {
      if (alpha_div_lambda * (log(K) + log((query_length-l) *
                             (db_length-num_of_db_seqs*l))) + beta >= l)
        length_adjustment = (GtUword) l;
    }
  }

  return length_adjustment;
}

GtUword gt_evalue_calculate_searchspace(const GtQuerymatch *querymatch,
                                        const GtEncseq *encseq,
                                        double alpha_div_lambda,
                                        double beta,
                                        double K)
{
  GtUword total_length_of_db,
          num_of_db_seqs,
          actual_db_length,
          actual_query_length,
          effective_db_length,
          effective_query_length,
          length_adjustment;

  total_length_of_db = gt_encseq_total_length(encseq);
  num_of_db_seqs = gt_encseq_num_of_sequences(encseq);
  /* db length without seperators */
  actual_db_length = total_length_of_db - num_of_db_seqs;

  actual_query_length = gt_querymatch_querylen(querymatch);

  length_adjustment = gt_evalue_calculate_length_adjustment(actual_query_length,
                                                            actual_db_length,
                                                            num_of_db_seqs,
                                                            alpha_div_lambda,
                                                            beta,
                                                            K);

  effective_query_length = actual_query_length - length_adjustment;
  effective_db_length = actual_db_length -
                                      (num_of_db_seqs * length_adjustment);

  return gt_safe_mult_ulong(effective_query_length, effective_db_length);
}

double gt_evalue_calculate(const GtKarlinAltschulStat *ka,
                           GtUword searchspace,
                           double bit_score)
{
  double logK, lambda;
  gt_assert(ka);

  lambda = gt_karlin_altschul_stat_get_lambda(ka);
  logK = gt_karlin_altschul_stat_get_logK(ka);

  return searchspace * exp(-lambda * bit_score + logK);
}
