/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/
#ifndef BIOFMI2_MISC_H_INCLUDED
#define BIOFMI2_MISC_H_INCLUDED

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif /* HAVE_STDINT_H */
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif /* HAVE_INTTYPES_H */
#ifndef DEBUG
#define NDEBUG
#endif /* DEBUG */
#include <assert.h>
#include <stdlib.h>

#include <libgtcore/minmax.h>

/**
 * \file biofmi2misc.h
 * \brief Simple helper routines for distribution computations.
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
  while(n > 1)
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
  if(k == 0 || k == n)
    return 1;
  else if(k < n/2)
    return binomialCoeff(n, n - k);
  else
  {
    unsigned long i = k;
    accum = ++i;
    while(i < n)
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
multinomialCoeff(uint32_t n, size_t numBins, const uint32_t binSizes[])
{
  unsigned long accum = 1, nfac;
  size_t i, maxBin = 0, maxBinSize = 0;
#ifndef NDEBUG
  unsigned long binSum = 0;
#endif
  assert(n > 0 && numBins > 0 && binSizes);
  for(i = 0; i < numBins; ++i)
  {
#ifndef NDEBUG
    binSum += binSizes[i];
#endif
    if(binSizes[i] > maxBinSize)
    {
      maxBinSize = binSizes[i];
      maxBin = i;
    }
  }
  assert(binSum <= n);
  for(nfac = maxBinSize + 1; nfac <= n; ++nfac)
    accum *= nfac;
  for(i = 0; i < numBins; ++i)
    if(i != maxBin)
      accum /= factorial(binSizes[i]);
  return accum;
}

/**
 * Compute number of digits a value would require when displayed in
 * selected number system (i.e. 2=binary, 10=decimal).
 * @param d value to be displayed
 * @param base number of symbols in output alphabet
 * @return number of digits required for output
 */
static inline int
digitPlaces(long d, int base)
{
  int l=1;
  while(d/=base)
    ++l;
  return l;
}

/*
 * functionality to layout data at correct alignment (=> fewer
 * individual mallocs)
 */
enum {
  MAX_ALIGN_REQUIREMENT = 8,
  MIN_ALIGN_REQUREMENT = 4,
};

static inline size_t
offsetAlign(size_t offset, size_t sizeOfVal2Align)
{
  size_t alignBase = MAX(MIN_ALIGN_REQUREMENT,
                         MIN(sizeOfVal2Align, MAX_ALIGN_REQUIREMENT));
  return (offset/alignBase + ((offset%alignBase)?1:0))*alignBase;
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


#endif /* BIOFMI2_MISC_H_INCLUDED */
