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

#include <assert.h>

#include "libgtcore/minmax.h"

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
factorial(int n)
{
  unsigned long k = 1;
  while (n > 1)
    k *= n--;
  return k;
}

/**
 * \brief Compute binomial coefficient \f[{n\choose k} = \frac{n!}{k!\cdot (n
 * - k)!}\f]
 * @param n
 * @param k
 * @return \f$n\choose{k}\f$
 */
static inline unsigned long
binomialCoeff(unsigned long n, unsigned long k)
{
  unsigned long accum;
  assert(k <= n);
  if (k == 0 || k == n)
    return 1;
  else if (k < n/2)
    return binomialCoeff(n, n - k);
  else
  {
    unsigned long i = k;
    accum = ++i;
    while (i < n)
      accum *= ++i;
    return accum /= factorial(n - k);
  }
}

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
static inline unsigned long
multinomialCoeff(unsigned n, size_t numBins, const unsigned binSizes[])
{
  unsigned long accum = 1, nfac;
  size_t i, maxBin = 0, maxBinSize = 0;
#ifndef NDEBUG
  unsigned long binSum = 0;
#endif
  assert(n > 0 && numBins > 0 && binSizes);
  for (i = 0; i < numBins; ++i)
  {
#ifndef NDEBUG
    binSum += binSizes[i];
#endif
    if (binSizes[i] > maxBinSize)
    {
      maxBinSize = binSizes[i];
      maxBin = i;
    }
  }
  assert(binSum <= n);
  for (nfac = maxBinSize + 1; nfac <= n; ++nfac)
    accum *= nfac;
  for (i = 0; i < numBins; ++i)
    if (i != maxBin)
      accum /= factorial(binSizes[i]);
  return accum;
}

static inline unsigned long long
iPow(unsigned long long x, unsigned i)
{
   unsigned long long result = 1;
   while (i)
   {
     if (i & 1)
       result *= x;
     x *= x;
     i >>= 1;
   }
   return result;
}

#endif
