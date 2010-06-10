/*
  Copyright (c) 2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
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

#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/log_api.h"
#include "match/eis-bwtseq-siop.h"
#include "match/eis-bwtseq.h"
#include "match/eis-encidxseq.h"
#include "match/eis-mrangealphabet-siop.h"
#include "match/eis-mrangealphabet.h"
#include "match/sarr-def.h"

#include "match/shu-match.h"

static double pmax(double M, /* M value should be explored by simulation ??? */
            unsigned long x,
            double p,
            unsigned long subjectLength,
            int *thresholdReached,
            double *ln_n_fac,
            double *s1,
            unsigned int n_s)
{

  unsigned long k;
  double s = 0, ln_x_choose_k;
  double ln, ln1, m1, m, delta;
  double pS = p;

  if (x > n_s)
  {
    /* change this to standard GT-behaviour XXX*/
    printf("Error: x = %lu."
           " The maximum number of elements in the array"
           " s1 should be increased!\n", x);
  }
  gt_assert(x <= n_s);
  if (s1[x] != 0)
  {
    return s1[x];
  }

  for (k = 0; k <= x; k++)
  {
    if (x == k)
      ln_x_choose_k = 0.0;
    else
      ln_x_choose_k =  ln_n_fac[x] - ln_n_fac[k] - ln_n_fac[x-k];
    gt_log_log("x = %lu, k = %lu", x, k);
    gt_log_log("ln_x_choose_k = %8.32f", ln_x_choose_k);

    m = (pow (2.0, (double) x) *
         pow (p, (double) k) *
         pow (0.5 - p, (double) x - k) *
         pow (1.0 - pow (pS, (double) k) * pow (0.5 - pS, (double) x - k),
              (double) subjectLength));
    /*gt_log_log("m: %f", m);*/
    if (m == 0)
    {
      delta = 0;
    } else if (m >= M)
    {
      ln = log(m);
      if (ln == -HUGE_VAL)
      {
        gt_log_log("huge!");
        delta = 0;
      } else
        delta = exp (ln + ln_x_choose_k);
    } else
    {
      gt_log_log("m kleiner als M");
      m1 = 1 + m;  /* for small values of m - to avoid overflow (-INF) */
      ln1 = log(m1);
      delta = exp (ln1 + ln_x_choose_k) - exp (ln_x_choose_k);
    }
    s += delta;
    if (s >= 1.0)
    {
      s = 1;
      *thresholdReached = 1;
      break;
    }
  }  /* end for */
  s1[x] = s;
  /*gt_log_log("pmax: %f", s);*/
  return s;
}

static double expShulen(double T, /* absolute error */
                 double M,
                 double d,
                 double p,
                 unsigned long subjectLength,
                 double *ln_n_fac,
                 double *s1,
                 unsigned int n_s)
{
  unsigned long i;
  int thresholdReached = 0;

  double prob_i, probOld, delta, factor;

  double e = 0.0;    /* expectation */
  double t = 1.0 - d;
  double p_t = t;

  probOld  = 0;

  /*since for i = 0, the whole expression is 0*/
  for (i = 1; i < subjectLength; i++)
  {
    factor = 1.0 - p_t;
    if (!thresholdReached)
    {
      prob_i = factor * pmax (M,
                              i,
                              p,
                              subjectLength,
                              &thresholdReached,
                              ln_n_fac,
                              s1,
                              n_s);
    } else
    {
      prob_i = factor;  /* prob_i = factor * s, where s = 1 */
    }
    delta = (prob_i - probOld) * i;  /* delta should always be positive */
    e += delta;    /* expectation of avg shulen length(Q, S) */
    /* check error */
    if (e >= 1 && delta / e <= T)
    {
      break;
    }
    p_t *= t;
    probOld = prob_i;
  }
  /*gt_log_log("expected shulength: %f", e);*/
  return e;
}

unsigned int gt_pck_getShuStringLength(const BWTSeq *bwtSubject,
                                       const GtUchar *query,
                                       size_t queryLength)
{
  const GtUchar *qptr, *qend;
  Symbol curChar;
  const MRAEnc *alphabet;
  unsigned long start, end;
  size_t retval;

  gt_assert(bwtSubject && query);
  alphabet = BWTSeqGetAlphabet(bwtSubject);

  qptr = query;
  qend = query + queryLength;

  curChar = MRAEncMapSymbol(alphabet, *qptr);
  /*gt_log_log("query[%lu]=%d",(unsigned long) (qptr-query),(int) *qptr);*/

  qptr++;
  start = bwtSubject->count[curChar];
  end = bwtSubject->count[curChar + 1];
  /*gt_log_log("start=%lu, end=%lu", start, end);*/
  for (/* Nothing */; start < end && qptr < qend; qptr++)
  {
    /*gt_log_log("query[%lu]=%d",(unsigned long) (qptr-query),(int) *qptr);*/
    GtUlongPair occPair;
    curChar = MRAEncMapSymbol(alphabet, *qptr);
    occPair = BWTSeqTransformedPosPairOcc(bwtSubject,
                                          curChar,
                                          start,
                                          end);
    start = bwtSubject->count[curChar] + occPair.a;
    end = bwtSubject->count[curChar] + occPair.b;
    /*gt_log_log("start=%lu, end=%lu", start, end);*/
  }
  if (qptr == qend && start < end)
    retval = queryLength + 1;
  else
    retval = (size_t) (qptr - query);
  return (unsigned int) retval;
}

double gt_pck_getGCcontent(const BWTSeq *bwtSubject,
                           const GtAlphabet *alphabet)
{
  unsigned long  c, length;
  double gc;
  const MRAEnc *FM_alphabet;
  GtUchar c_sym;

  FM_alphabet = BWTSeqGetAlphabet(bwtSubject);

  c_sym = MRAEncMapSymbol(FM_alphabet,
                          gt_alphabet_encode(alphabet, 'c'));
  length = bwtSubject->seqIdx->seqLen;

  c = bwtSubject->count[c_sym+1] - bwtSubject->count[c_sym];

  gc = c * 2 / (double) (length - 2);
  return gc;
}

/* calculate divergence */
double gt_divergence(double E, /* relative error for shulen length */
                   double T, /* absolute error */
                   double M,
                   double shulen,
                   unsigned long subjectLength,
                   double gc,
                   double *ln_n_fac,
                   unsigned int n_s)
{
  double p, q;
  double du, dl, dm, t, d;
  double *s1;
  s1 = gt_calloc(n_s + 1, sizeof (double));

  p = gc / 2;
  q = (1.0 - gc) / 2.0;
  du = 0;
  dl = 1.0 - (2 * p * p + 2 * q * q);  /* dl < 0.75 */
  /*this should become user definable*/
  t = THRESHOLD;

  while ((dl - du) / 2.0 > t)
  {
    dm = (du + dl) / 2.0;
    if (shulen < expShulen (T, M, dm, p, subjectLength, ln_n_fac, s1, n_s))
    {
      du = dm;
    } else
    {
      dl = dm;
    }
    /* test the relative error between du and dl; if it is smaller than some
     * threshold, then break !! */
    if (fabs (dl - du) / dl <= E)
    {
      break;
    }
  }
  d = (du + dl) / 2.0;
  /*gt_log_log("divergence: %2f", d);*/
  gt_free(s1);
  return d;
}

double *gt_get_ln_n_fac(int n)
{
  int i;
  double *ln_n_fac;

  ln_n_fac = gt_calloc(n + 1, sizeof (double));
  gt_assert(ln_n_fac != NULL);

  for (i = 0; i <= n; i++)
  {
    if (i == 0)
      ln_n_fac[i] = log(1);
    else
      ln_n_fac[i] = log(i) + ln_n_fac[i-1];
  }
  return ln_n_fac;
}

double gt_calculateKr(double d)
{
  double kr;

  kr = -0.75 * log (1 - 4.0 / 3.0 * d);

  return kr;
}
