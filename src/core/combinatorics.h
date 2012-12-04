/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include "core/assert_api.h"
#include "core/error_api.h"
#include "core/minmax.h"
#include "core/unused_api.h"

#ifdef _LP64
#define GT_BINOMIAL_MAX_N_LN 66UL
#define GT_BINOMIAL_MAX_N 62UL
#define GT_BINOMIAL_MAX_N_DP 67UL
#else
#define GT_BINOMIAL_MAX_N_LN 32UL
#define GT_BINOMIAL_MAX_N 30UL
#define GT_BINOMIAL_MAX_N_DP 33UL
#endif

/**
 * \file combinatorics.h
 * \brief Simple routines for distribution computations relevant to
 * combinatorial problems.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

/**
 * Computes \f$n! = 1 \cdot 2 \cdot \dots \cdot (n - 1) \cdot n\f$
 * @param n number of which to compute factorial
 * @return \f$n!\f$
 */
static inline unsigned long
factorial(unsigned n)
{
  unsigned long k = 1UL;
  while (n > 1U)
    k *= n--;
  return k;
}

/* calculate n choose k using exp(ln(n!) - ln(k!) - ln((n-k)!)). Returned value
   might deviate from correct result for large n. Overflows for
   n > GT_BINOMIAL_MAX_N_LN */
unsigned long gt_binomialCoeff_with_ln(unsigned long n, unsigned long k);

/* calculate n choose k using a dp table. Overflows for
   n> GT_BINOMIAL_MAX_N_DP */
unsigned long gt_binomialCoeff_dp(unsigned long n, unsigned long k);

/* naive implementation of n choose k, but already somewhat optimised.
   Overflows for n> GT_BINOMIAL_MAX_N */
unsigned long gt_binomialCoeff(unsigned long n, unsigned long k);

/**
 * \brief Compute multinomial coefficient
 * \f[{n\choose k_1, k_2,\dots,k_m} = \frac{n!}{k_1!\cdot
 * k_2!\cdot\dots\cdot k_m!}\f] where \f$m\f$ equals \link numBins\endlink
 * @param n
 * @param numBins
 * @param binSizes points to array containing \link numBins\endlink values which
 * represent the \f$k_i\f$
 * @return \f$n\choose{k_1, k_2,\dots,k_m}\f$
 */
/*@unused@*/ static inline unsigned long
multinomialCoeff(unsigned n, size_t numBins, const unsigned binSizes[])
{
  unsigned long accum = 1UL, nfac;
  size_t i, maxBin = 0, maxBinSize = 0;
#ifndef NDEBUG
  unsigned long binSum = 0;
#endif
  gt_assert(n > 0 && numBins > 0 && binSizes);
  for (i = 0; i < numBins; ++i) {
#ifndef NDEBUG
    binSum += binSizes[i];
#endif
    if ((size_t) binSizes[i] > maxBinSize) {
      maxBinSize = (size_t) binSizes[i];
      maxBin = i;
    }
  }
  gt_assert(binSum <= (unsigned long) n);
  for (nfac = (unsigned long) maxBinSize + 1; nfac <= (unsigned long) n; ++nfac)
    accum *= nfac;
  for (i = 0; i < numBins; ++i)
    if (i != maxBin)
      accum /= factorial(binSizes[i]);
  return accum;
}

/*@unused@*/ static inline unsigned long long
iPow(unsigned long long x, unsigned i)
{
   unsigned long long result = 1ULL;
   while (i) {
     if (i & 1)
       result *= x;
     x *= x;
     i >>= 1;
   }
   return result;
}

int gt_combinatorics_unit_test(GT_UNUSED GtError *err);
void gt_combinatorics_init(void);
void gt_combinatorics_clean(void);
#endif
