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
#include "core/safearith.h"

#ifdef _LP64
#define GT_BINOMIAL_MAX_N_LN 66UL
#define GT_BINOMIAL_MAX_N 62UL
#define GT_BINOMIAL_MAX_N_DP 67UL
#else
#define GT_BINOMIAL_MAX_N_LN 32UL
#define GT_BINOMIAL_MAX_N 30UL
#define GT_BINOMIAL_MAX_N_DP 33UL
#endif

/* Module Combinatorics includes some standard aproaches for combinatorical
   problems. Some functions have efficient implementations.
   Some of the simpler functions are meant to be used inline. */

/* Returns n! = 1 * 2 * ... * (n - 1) * n
   <n> number of which to compute factorial */
static inline unsigned long gt_combinatorics_factorial(unsigned n);

/* calculate n choose k using exp(ln(n!) - ln(k!) - ln((n-k)!)). Returned value
   might deviate from correct result for large n. Overflows for
   n > GT_BINOMIAL_MAX_N_LN */
unsigned long                    gt_combinatorics_binomial_ln(unsigned long n,
                                                              unsigned long k);

/* calculate n choose k using a dp table. Overflows for
   n> GT_BINOMIAL_MAX_N_DP */
unsigned long                    gt_combinatorics_binomial_dp(unsigned long n,
                                                              unsigned long k);

/* naive implementation of n choose k, but already somewhat optimised.
   Overflows for n> GT_BINOMIAL_MAX_N */
unsigned long                    gt_combinatorics_binomial_simple(
                                                               unsigned long n,
                                                               unsigned long k);

/* Returns multinomial coefficient
   n choose k_1, k_2 ... k_m = n! / k_1! * k_2! * ... * k_m!
   <numBins> equals m, <binSizes> points to array which contains the k_i. */
static inline unsigned long      gt_combinatorics_multinomial(
                                                     unsigned n,
                                                     size_t numBins,
                                                     const unsigned binSizes[]);

/* Returns <x>^<i> */
static inline unsigned long long gt_combinatorics_i_pow(unsigned long long x,
                                                        unsigned i);

int                              gt_combinatorics_unit_test(
                                                        GT_UNUSED GtError *err);
void                             gt_combinatorics_init(void);
void                             gt_combinatorics_clean(void);

#include "core/combinatorics_impl.h"
#endif
