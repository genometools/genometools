/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include <float.h>

#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/log_api.h"
#include "core/mathsupport.h"
#include "match/shu-divergence.h"

static double pmax(double M, /* M value should be explored by simulation ??? */
            unsigned long x,
            double p,
            unsigned long subjectLength,
            int *thresholdReached,
            double *ln_n_fac,
            double *s1,
            unsigned long n_s)
{

  unsigned long k;
  double s = 0.0, ln_x_choose_k;
  double ln, ln1, m1, m, delta;

  if (x > n_s) {
    /* change this to standard GT-behaviour XXX*/
    printf("Error: x = %lu."
           " The maximum number of elements in the array"
           " s1 should be increased!\n", x);
  }
  gt_assert(x <= n_s);
  if (s1[x] != 0.0)
    return s1[x];

  for (k = 0; k <= x; k++) {
    double m_a, m_b, m_c, m_d, m_e;
    if (x == k)
      ln_x_choose_k = 0.0;
    else
      ln_x_choose_k =  ln_n_fac[x] - ln_n_fac[k] - ln_n_fac[x-k];

    m_a = pow (2.0, (double) x);
    m_b = pow (p, (double) k);
    m_c = pow (0.5 - p, (double) x - k);
    m_d = pow (p, (double) k);
    m_e = pow (0.5 - p, (double) x -k);
    m = (m_a * m_b * m_c * pow (1.0 - m_d * m_e, (double) subjectLength));
    /* this is ok even with double, because of next if!*/
    if (m == 0.0)
      delta = 0.0;
    else if (M <= m) {
      ln = log(m);
      if (ln == -HUGE_VAL)
        delta = 0.0;
      else
        delta = exp (ln + ln_x_choose_k);
    }
    else {
      double delta_a, delta_b;
      m1 = 1 + m;  /* for small values of m - to avoid overflow (-INF) */
      ln1 = log(m1);
      delta_a = exp (ln1 + ln_x_choose_k);
      delta_b = exp (ln_x_choose_k);
      delta = delta_a - delta_b;
    }
    s += delta;
    if (s >= 1.0) {
      s = 1.0;
      *thresholdReached = 1;
      break;
    }
  }  /* end for */
  s1[x] = s;
  return s;
}

static double expShulen(double T, /* absolute error */
                 double M, /* min logarithm */
                 double d,
                 double p,
                 unsigned long subjectLength,
                 double *ln_n_fac,
                 double *s1,
                 unsigned long n_s)
{
  unsigned long i;
  int thresholdReached = 0;

  double prob_i, probOld, delta, factor;

  double e = 0.0;    /* expectation */
  double t = 1.0 - d;
  double p_t = t;

  probOld  = 0.0;

  /*since for i = 0, the whole expression is 0*/
  for (i = 1LU; i < subjectLength; i++) {
    factor = 1.0 - p_t;
    if (!thresholdReached) {
      prob_i = factor * pmax(M, i, p,
                             subjectLength,
                             &thresholdReached,
                             ln_n_fac,
                             s1,
                             n_s);
    }
    else
      prob_i = factor;  /* prob_i = factor * s, where s = 1 */
    delta = (prob_i - probOld) * i;  /* delta should always be positive */
    gt_assert(gt_double_equals_double(delta, 0.0) ||
              gt_double_larger_double(delta, 0.0));
    e += delta;    /* expectation of avg shulen length(Q, S) */
    /* check error */
    if (1.0 <= e && delta / e <= T)
      break;
    p_t *= t;
    probOld = prob_i;
  }
  return e;
}

/* calculate divergence */
double gt_divergence(double E, /* relative error for divergence calculation */
                   double T, /* absolute error for exp shulen*/
                   double M, /* minimum for logarithm */
                   double threshold, /* abs error for divergence */
                   double shulen,
                   unsigned long subjectLength,
                   double gc,
                   double *ln_n_fac,
                   unsigned long n_s)
{
  double p, q,
         du, dl, dm, d, exp_shulen;
  double *s1;
  s1 = gt_calloc((size_t) n_s + 1, sizeof (double));

  p = gc / 2;
  q = (1.0 - gc) / 2.0;
  du = 0.0;
  dl = 1.0 - (2 * p * p + 2 * q * q);  /* dl < 0.75 */

  while (gt_double_smaller_double(threshold, (dl - du) / 2.0)) {
    dm = (du + dl) / 2.0;
    exp_shulen = expShulen(T, M, dm, p, subjectLength, ln_n_fac, s1, n_s);
    if (gt_double_smaller_double(shulen, exp_shulen))
      du = dm;
    else
      dl = dm;
    /* test the relative error between du and dl; if it is smaller than some
     * threshold, then break !! */
    if (fabs (dl - du) / dl <= E)
      break;
  }
  d = (du + dl) / 2.0;
  gt_free(s1);
  return d;
}

double *gt_get_ln_n_fac(unsigned long n)
{
  unsigned long i;
  double *ln_n_fac;

  ln_n_fac = gt_calloc((size_t) n + 1, sizeof (double));

  ln_n_fac[0] = log(1.0);
  for (i = 1UL; i <= n; i++) {
    ln_n_fac[i] = log((double) i) + ln_n_fac[i-1];
  }
  return ln_n_fac;
}

double gt_calculateKr(double d)
{
  double kr = -0.75 * log (1 - 4.0 / 3.0 * d);
  return kr;
}
