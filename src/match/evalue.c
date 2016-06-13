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
