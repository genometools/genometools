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

#include <float.h>
#include <limits.h>
#include <math.h>

#include "core/array2dim_api.h"
#include "core/combinatorics.h"
#include "core/divmodmul.h"
#include "core/ensure.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/safearith.h"
#include "core/warning_api.h"

#define UNDEFTABVALUE ULONG_MAX

static double *ln_n_fac_tab = NULL;
static unsigned long **binomial_dp_tab = NULL;

static inline void init_ln_n_fac_tab(void)
{
  unsigned long idx;

  if (ln_n_fac_tab == NULL) {
    ln_n_fac_tab = gt_calloc((size_t) (GT_BINOMIAL_MAX_N_LN + 1),
                             sizeof (double));

    ln_n_fac_tab[0] = log(1.0);
    for (idx = 1UL; idx <= GT_BINOMIAL_MAX_N_LN; idx++) {
      ln_n_fac_tab[idx] = log((double) idx) + ln_n_fac_tab[idx-1];
    }
    gt_log_log("ln_fac_max: %lf", ln_n_fac_tab[GT_BINOMIAL_MAX_N_LN]);
  }
}

static inline void init_binomial_dp_tab(void)
{
  if (binomial_dp_tab == NULL) {
    unsigned long idx, jdx,
                  rows = GT_BINOMIAL_MAX_N_DP + 1UL,
                  cols = GT_DIV2(GT_BINOMIAL_MAX_N_DP) + 1;
    gt_array2dim_malloc(binomial_dp_tab, rows, cols);
    for (idx=0; idx < rows; idx++)
      for (jdx=0; jdx < cols; jdx++)
        binomial_dp_tab[idx][jdx] = UNDEFTABVALUE;
  }
}

void gt_combinatorics_init(void)
{
  init_ln_n_fac_tab();
  init_binomial_dp_tab();
}

void gt_combinatorics_clean(void)
{
  gt_free(ln_n_fac_tab);
  ln_n_fac_tab = NULL;
  gt_array2dim_delete(binomial_dp_tab);
  binomial_dp_tab = NULL;
}

static inline unsigned long gt_combinatorics_binomial_dp_rec(unsigned long n,
                                                             unsigned long k)
{
  if (k == 0 || n <= k)
    return 1UL;
  if (n == 0)
    return 0;
  if (binomial_dp_tab[n][k] == UNDEFTABVALUE) {
    if (binomial_dp_tab[n - 1][k - 1] == UNDEFTABVALUE)
      binomial_dp_tab[n - 1][k - 1] =
        gt_combinatorics_binomial_dp_rec(n - 1, k - 1);
    if (binomial_dp_tab[n - 1][k] == UNDEFTABVALUE)
      binomial_dp_tab[n - 1][k] = gt_combinatorics_binomial_dp_rec(n - 1, k);
    gt_safe_add(binomial_dp_tab[n][k],
                binomial_dp_tab[n - 1][k - 1],
                binomial_dp_tab[n - 1][k]);
  }
  return binomial_dp_tab[n][k];
}

unsigned long gt_combinatorics_binomial_dp(unsigned long n, unsigned long k)
{
  if (k == 0 || n <= k)
    return 1UL;
  if (n == 0)
    return 0;
  gt_assert(n <= GT_BINOMIAL_MAX_N_DP);
  if (GT_DIV2(n) < k)
    k = n - k;
  return gt_combinatorics_binomial_dp_rec(n, k);
}

unsigned long gt_combinatorics_binomial_ln(unsigned long n, unsigned long k)
{
  if (k == 0 || n <= k)
    return 1UL;
  if (n < k)
    return 0;
  gt_assert(n <= GT_BINOMIAL_MAX_N_LN);
  if (GT_DIV2(n) < k)
    k = n - k;

  return gt_safe_cast2ulong(gt_round_to_long(exp(ln_n_fac_tab[n] -
                                                 ln_n_fac_tab[k] -
                                                 ln_n_fac_tab[n - k])));
}

unsigned long gt_combinatorics_binomial_simple(unsigned long n, unsigned long k)
{
  unsigned long idx;
  unsigned long result;
  if (n < k)
    return 0;
  if (k == 0 || k == n)
    return 1UL;
  gt_assert(n <= GT_BINOMIAL_MAX_N);
  if (GT_DIV2(n) < k)
    k = n - k;
  result = n - k + 1;
  for (idx = 2UL; idx <= k; idx++) {
#ifdef _LP64
    result = (unsigned long) gt_safe_mult_u64((uint64_t) result,
                                              (uint64_t) (n - k + idx));
#else
    result = (unsigned long) gt_safe_mult_u32((uint32_t) result,
                                              (uint32_t) (n - k + idx));
#endif
    result /= idx;
  }
  return result;
}

int gt_combinatorics_unit_test(GtError *err)
{
  int had_err = 0;
  unsigned long n,k;
  static const unsigned long max_n = GT_BINOMIAL_MAX_N,
                             max_fac_stable = 47UL;
  gt_error_check(err);
  for (n = 0; n <= max_n; n++) {
    for (k = 0; k <= GT_DIV2(n); k++) {
      unsigned long a = gt_combinatorics_binomial_dp(n, k),
                    b = gt_combinatorics_binomial_simple(n, k),
                    c;
      gt_ensure(had_err, a == b);
      if (n < max_fac_stable) {
        c = gt_combinatorics_binomial_ln(n,k);
        gt_ensure(had_err, c == a);
      }
    }
  }
  return had_err;
}
