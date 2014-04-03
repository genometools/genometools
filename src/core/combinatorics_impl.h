/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef COMBINATORICS_IMPL_H
#define COMBINATORICS_IMPL_H

#include "core/assert_api.h"
#include "core/error_api.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/safearith.h"

#if defined (_LP64) || defined (_WIN64)
#define GT_BINOMIAL_MAX_N_LN 66UL
#define GT_BINOMIAL_MAX_N 62UL
#define GT_BINOMIAL_MAX_N_DP 67UL
#else
#define GT_BINOMIAL_MAX_N_LN 32UL
#define GT_BINOMIAL_MAX_N 30UL
#define GT_BINOMIAL_MAX_N_DP 33UL
#endif

static inline GtUword gt_combinatorics_factorial(unsigned n)
{
  GtUword k = 1UL;
  while (n > 1U)
    k *= n--;
  return k;
}

/*@unused@*/
static inline GtUint64 gt_combinatorics_i_pow(GtUint64 x,
                                              unsigned int i)
{
   GtUint64 result = 1ULL;
   if (x == 2ULL)
     return 1ULL << (GtUint64) i;
   while (i != 0) {
     if (i & 1)
       result *= x;
     x *= x;
     i >>= 1;
   }
   return result;
}

/*@unused@*/
static inline GtUword gt_combinatorics_multinomial(unsigned n,
                                                   size_t numBins,
                                                   const unsigned binSizes[])
{
  GtUword accum = 1UL, nfac;
  size_t i, maxBin = 0, maxBinSize = 0;
#ifndef NDEBUG
  GtUword binSum = 0;
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
  gt_assert(binSum <= (GtUword) n);
  for (nfac = (GtUword) maxBinSize + 1; nfac <= (GtUword) n; ++nfac)
    accum *= nfac;
  for (i = 0; i < numBins; ++i)
    if (i != maxBin)
      accum /= gt_combinatorics_factorial(binSizes[i]);
  return accum;
}

#endif
