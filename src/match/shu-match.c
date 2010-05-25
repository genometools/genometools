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
#include <gsl/gsl_sf_gamma.h>
#include "gsl/gsl_nan.h"

#include "core/assert_api.h"
#include "core/log_api.h"
#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-siop.h"
#include "match/eis-encidxseq.h"
#include "match/eis-mrangealphabet-siop.h"
#include "match/eis-mrangealphabet.h"
#include "match/sarr-def.h"

#include "match/shu-match.h"

double pmax(double M,
            unsigned long x,
            double p,
            unsigned long queryLength,
            int *thresholdReached,
            double **lnChoose,
            int n_choose,
            double *s1);

double expShulen(double T,
                 double M,
                 double d,
                 double p,
                 unsigned long queryLength,
                 double **lnChoose,
                 int n_choose,
                 double *s1);

unsigned int gt_pck_getShuStringLength(const BWTSeq *bwtSubject,
                                       const GtUchar *query,
                                       size_t queryLength)
{
  const GtUchar *qptr, *qend;
  GtUchar curChar;
  const MRAEnc *alphabet;
  unsigned long start, end;
  size_t retval;

  gt_assert(bwtSubject && query);
  alphabet = BWTSeqGetAlphabet(bwtSubject);

  qptr = query;
  qend = query + queryLength;

  curChar = MRAEncMapSymbol(alphabet, *qptr);
  gt_log_log("query[%lu]=%d",(unsigned long) (qptr-query),(int) *qptr);

  qptr += 1;
  start = bwtSubject->count[curChar];
  end = bwtSubject->count[curChar + 1];
  gt_log_log("start=%lu, end=%lu", start, end);
  while (start < end && qptr != qend)
  {
    gt_log_log("query[%lu]=%d",(unsigned long) (qptr-query),(int) *qptr);
    GtUlongPair occPair;
    curChar = MRAEncMapSymbol(alphabet, *qptr);
    occPair = BWTSeqTransformedPosPairOcc(bwtSubject,
                                          (Symbol) curChar,
                                          start,
                                          end);
    start = bwtSubject->count[curChar] + occPair.a;
    end = bwtSubject->count[curChar] + occPair.b;
    gt_log_log("start=%lu, end=%lu", start, end);
    qptr += 1;
  }
  if (qptr == qend && start < end)
    retval = queryLength + 1;
  else
    retval = qptr - query;
  return (unsigned int) retval;
}

double gt_pck_getGCcontent(const BWTSeq *bwtSubject,
                           const GtAlphabet *alphabet)
{
  unsigned long /*a,*/ c,/* g, t,*/ length;
  double gc;
  const MRAEnc *FM_alphabet;
  GtUchar /*a_sym,*/ c_sym/*, g_sym, t_sym*/;

  FM_alphabet = BWTSeqGetAlphabet(bwtSubject);

/*  a_sym = MRAEncMapSymbol(FM_alphabet,
                          gt_alphabet_encode(alphabet, 'a')); */
  c_sym = MRAEncMapSymbol(FM_alphabet,
                          gt_alphabet_encode(alphabet, 'c'));
/*  g_sym = MRAEncMapSymbol(FM_alphabet,
                          gt_alphabet_encode(alphabet, 'g')); */
/*  t_sym = MRAEncMapSymbol(FM_alphabet,
                          gt_alphabet_encode(alphabet, 't')); */
  length = bwtSubject->seqIdx->seqLen;

/*  a = bwtSubject->count[a_sym+1] - bwtSubject->count[a_sym]; */
  c = bwtSubject->count[c_sym+1] - bwtSubject->count[c_sym];
/*  g = bwtSubject->count[g_sym+1] - bwtSubject->count[g_sym]; */
/*  t = bwtSubject->count[t_sym+1] - bwtSubject->count[t_sym]; */

/*  gt_log_log("acgt: %lu %lu %lu %lu",
             a, c, g, t); */
  gt_log_log("%lu / %lu",
             c, length);
  /*double stranded: g = c */
  gc = c * 2 / (double) length;
  return gc;
}

/* calculate divergence */
double divergence (double E, /* relative error for shulen length */
                   double T, /* absolute error */
                   double M,
                   double shulen,
                   unsigned long queryLength,
                   double gc,
                   int n_choose, /* change that in future
                                  * to be user spezified */
                   double **lnChoose)
{
  double p, q;
  double du, dl, dm, t, d;
  double s1[N_S + 1] = { 0 };

  p = gc / 2;
  q = (1.0 - gc) / 2.0;
  du = 0;
  dl = 1.0 - (2 * p * p + 2 * q * q);  /* dl < 0.75 */
  /*this should become user definable*/
  t = THRESHOLD;

  while ((dl - du) / 2.0 > t)
  {
    dm = (du + dl) / 2.0;
    if (shulen < expShulen (T, M, dm, p, queryLength, lnChoose, n_choose, s1))
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
  return d;
}

double expShulen(double T, /* absolute error */
                 double M,
                 double d,
                 double p,
                 unsigned long queryLength,
                 double **lnChoose,
                 int n_choose,
                 double *s1)
{
  unsigned long i;
  int thresholdReached = 0;

  double prob_i, probOld, delta, factor;

  double e = 0.0;    /* expectation */
  double t = 1.0 - d;
  double p_t = t;

  probOld  = 0;

  /*since for i = 0, the whole expression is 0*/
  for (i = 1; i < queryLength; i++)
  {
    factor = 1.0 - p_t;
    if (!thresholdReached)
    {
      prob_i = factor * pmax (M,
                              i,
                              p,
                              queryLength,
                              &thresholdReached,
                              lnChoose,
                              n_choose,
                              s1);
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
  return e;
}

double pmax(double M, /* M value should be explored by simulation ??? */
            unsigned long x,
            double p,
            unsigned long queryLength,
            int *thresholdReached,
            double **lnChoose,
            int n_choose,
            double *s1)
{

  unsigned long k;
  double s = 0, x_choose_k;
  double t, t1, m, delta;
  double pS = p;

  if (x > N_S)
  {
    /* change this to standard GT-behaviour XXX*/
    printf("Error: x = %lu."
           " The maximum number of elements in the array"
           " s1 should be increased!\n", x);
  }
  if (s1[x] != 0)
  {
    return s1[x];
  }

  s = 0;
  for (k = 0; k <= x; k++)
  {
    if (x < n_choose)
    {
      if (lnChoose[x][k] == 0)
      {
        lnChoose[x][k] = gsl_sf_lnchoose (x, k);
      }
      x_choose_k = lnChoose[x][k];
    } else
    {
      x_choose_k = gsl_sf_lnchoose (x, k);
    }
    m = (pow (2.0, (double) x) *
         pow (p, (double) k) *
         pow (0.5 - p, (double) x - k) *
         pow (1.0 - pow (pS, (double) k) * pow (0.5 - pS, (double) x - k),
              (double) queryLength));
    if (m == 0)
    {
      delta = 0;
    } else if (m >= M)
    {
      t = log (m);
      if (t == GSL_NEGINF)
      {
        delta = 0;
      } else
      {
        delta = exp (t + x_choose_k);
      }
    } else
    {
      t1 = log (1 + m);  /* for small values of m - to avoid overflow (-INF) */
      delta = exp (t1 + x_choose_k) - exp (x_choose_k);
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
  return s;
}

double **initializeLnChoose(unsigned long dim)
{
  double **as = NULL;    /* matrix N x N */
  unsigned long i, j;

  /* XXX gt_malloc benutzen*/
  as = malloc (dim * sizeof (double *));
  for (i = 0; i < dim; i++)
  {
    as[i] = malloc (dim * sizeof (double));
  }
  /* initialize ir elements */
  for (i = 0; i < dim; i++)
  {
    for (j = 0; j < dim; j++)
    {
      as[i][j] = 0.0;
    }
  }
  return as;
}

/* free matrix */
void freeLnChoose(double **matrix,
                  unsigned long dim)
{
  unsigned long i;
  for (i = 0; i < dim; i++)
    {
      free (matrix[i]);
    }
  free (matrix);
}

double calculateKr(double d)
{
  double kr;

  kr = -0.75 * log (1 - 4.0 / 3.0 * d);

  return kr;
}
