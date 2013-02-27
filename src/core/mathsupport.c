/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2006-2012 Center for Bioinformatics, University of Hamburg

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
#include <stdio.h>
#include <inttypes.h>
#include <limits.h>
#include "core/ensure.h"
#include "core/mathsupport.h"
#include "core/yarandom.h" /* necessary to define random() correctly */

#define GT_DBL_MAX_ABS_ERROR 1.0E-100
#define GT_DBL_MAX_REL_ERROR 1.0E-8

#ifndef S_SPLINT_S
double gt_logsum(double p1, double p2)
{
  if (p1 > p2)
    return (p1-p2 > 50.0) ? p1 : p1 + log(1.0 + exp(p2-p1));
  return (p2-p1 > 50.0) ? p2 : p2 + log(1.0 + exp(p1-p2));
}

static inline bool gt_double_relative_equal(double d1, double d2)
{
  double relerr;
  if (fabs(d1 - d2) < GT_DBL_MAX_ABS_ERROR)
    return true;
  if (fabs(d2) > fabs(d1))
    relerr = fabs((d1 - d2) / d2);
  else
    relerr = fabs((d1 - d2) / d1);
  if (relerr <= GT_DBL_MAX_REL_ERROR)
    return true;
  return false;
}

bool gt_double_equals_one(double d)
{
  return gt_double_relative_equal(d, 1.0);
}

bool gt_double_equals_double(double d1, double d2)
{
  return gt_double_relative_equal(d1, d2);
}

int gt_double_compare(double d1, double d2)
{
  if (gt_double_relative_equal(d1, d2))
    return 0;
  if (d1 > d2)
    return 1;
  return -1;
}

bool gt_double_smaller_double(double d1, double d2)
{
  return gt_double_compare(d1, d2) < 0;
}

bool gt_double_larger_double(double d1, double d2)
{
  return gt_double_compare(d1, d2) > 0;
}

unsigned long gt_rand_max(unsigned long maximal_value)
{
  unsigned long r;
  gt_assert(maximal_value);
  r = ((double) random() / ((double) RAND_MAX + 1) * (maximal_value + 1));
  gt_assert(r <= maximal_value);
  return r;
}

double gt_rand_max_double(double maximal_value)
{
  double r;
  gt_assert(maximal_value);
  r = ((double) random() / RAND_MAX) * maximal_value; /* XXX */
  gt_assert(r >= 0.0 && r <= maximal_value);
  return r;
}

double gt_rand_0_to_1(void)
{
  double r;
  r = (double) random() / RAND_MAX;
  gt_assert(r >= 0.0 && r <= 1.0);
  return r;
}

char gt_rand_char(void)
{
  int offset;
  char c;
  offset = gt_rand_max(25);
  c = 97 + offset;
  gt_assert(c >= 'a' && c <= 'z');
  return c;
}

/*
  Find the log base 2 of an integer in O(wordsize) operations.
  where N is the number of bits. There are faster methods, see
  \url{http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious}
*/
unsigned int gt_determinebitspervalue(unsigned long maxvalue)
{
  unsigned int bits = 0;
  unsigned long value;

  for (value = maxvalue; value > 0; value >>= 1)
  {
    bits++;
  }
#define GT_MAXLOG2VALUE 63
  gt_assert(bits <= GT_MAXLOG2VALUE);
#ifdef USEbuiltin_clzl
  {
    unsigned int fastbits = (unsigned int) (sizeof(unsigned long) * CHAR_BIT) -
                            __builtin_clzl(maxvalue);
    if (bits != fastbits)
    {
      fprintf(stderr,"maxvalue=%lu: bits = %u != %d = clzl\n",
                        maxvalue,bits,fastbits);
      exit(EXIT_FAILURE);
    }
  }
#endif
  return bits;
}

unsigned long gt_power_for_small_exponents(unsigned int base,
                                           unsigned int exponent)
{
  unsigned int logvalue = 0;

  if (exponent == 0) {
    return 1UL;
  }
  switch (base)
  {
    case 2UL:
      logvalue = 1U;
      break;
    case 4UL:
      logvalue = 2U;
      break;
    case 8UL:
      logvalue = 3U;
      break;
    case 16UL:
      logvalue = 4U;
      break;
    case 32UL:
      logvalue = 5U;
      break;
    case 64:
      logvalue = 6U;
      break;
    case 128:
      logvalue = 7U;
      break;
    case 256:
      logvalue = 8U;
      break;
  }
  if (logvalue > 0) {
    gt_assert(logvalue * exponent < sizeof (void *) * 8);
    return 1UL << (logvalue * exponent);
  } else {
    unsigned long powervalue = (unsigned long) base;

    while (exponent > 1U) {
      powervalue *= base;
      exponent--;
    }
    return powervalue;
  }
}

static inline double gt_round(double x)
{
  double integer = ceil(x);
  if (0.0 - x < GT_DBL_MAX_ABS_ERROR)
    return integer - x > 0.5 ? integer - 1.0 : integer;
  return integer - x >= 0.5 ? integer - 1.0 : integer;
}

long int gt_round_to_long(double x)
{
  int64_t intgr;
  double rounded = gt_round(x);
  intgr = (int64_t) rounded;
  /* If the fractional part is exactly 0.5, we need to check whether
  the rounded result is even. If it is not we need to add 1 to
  negative values and subtract one from positive values. */
  if ((fabs(((double) intgr) - x) == 0.5) & intgr) {
    /* 1 with the sign of result, i.e. -1 or 1. */
    intgr -= ((intgr >> 62) | ((int64_t) 1));
  }
  return (long int) intgr;
}

/* Make some unit tests for this? */
void gt_out_power_for_small_exponents(void)
{
  unsigned int exponent;

  for (exponent=1U; exponent<64U; exponent++) {
    printf("pow(2UL,%u)=%lu\n",exponent,
            gt_power_for_small_exponents(2U,exponent));
  }
  for (exponent=1U; exponent<32U; exponent++) {
    printf("pow(4UL,%u)=%lu\n",exponent,
            gt_power_for_small_exponents(4U,exponent));
  }
  for (exponent=1U; exponent<16U; exponent++) {
    printf("pow(8UL,%u)=%lu\n",exponent,
            gt_power_for_small_exponents(8U,exponent));
  }
  for (exponent=1U; exponent<32U; exponent++) {
    printf("pow(3UL,%u)=%lu\n",exponent,
            gt_power_for_small_exponents(3U,exponent));
  }
}

int gt_mathsupport_unit_test(GtError *err)
{
  int had_err = 0;
  double less_than_epsilon = 0.0000000000000001;
  gt_error_check(err);

  gt_ensure(had_err, !gt_double_equals_one(1.1));
  gt_ensure(had_err, gt_double_equals_one(1));
  gt_ensure(had_err, gt_double_equals_one(1+less_than_epsilon));
  gt_ensure(had_err, !gt_double_equals_one(-1-less_than_epsilon));

  gt_ensure(had_err, !gt_double_equals_double(1.0, 2.0));
  gt_ensure(had_err, !gt_double_equals_double(-1.0, 1.0));
  gt_ensure(had_err, !gt_double_equals_double(1.0, -1.0));
  gt_ensure(had_err, !gt_double_equals_double(-1.0, 1+less_than_epsilon));
  gt_ensure(had_err, !gt_double_equals_double(1.0, 1.1));
  gt_ensure(had_err, gt_double_equals_double(1.0, 1+less_than_epsilon));
  gt_ensure(had_err, gt_double_equals_double(1.0, 1.0));
  gt_ensure(had_err, gt_double_equals_double(0.0, 0.0));
  gt_ensure(had_err, gt_double_equals_double(-1.0, -1.0));
  gt_ensure(had_err, gt_double_equals_double(-1.0+less_than_epsilon, -1.0));
  gt_ensure(had_err, gt_double_equals_double(-1.0, -1.0+less_than_epsilon));
  gt_ensure(had_err, gt_double_equals_double(1.0+less_than_epsilon, 1.0));
  gt_ensure(had_err, gt_double_equals_double(1.0, 1.0+less_than_epsilon));

  gt_ensure(had_err, gt_double_compare(1.0, 1.0) == 0);
  gt_ensure(had_err, gt_double_compare(1.0, 1.1) < 0);
  gt_ensure(had_err, gt_double_compare(1.1, 1.0) > 0);
  gt_ensure(had_err, gt_double_compare(1.1, -1.0) > 0);
  gt_ensure(had_err, gt_double_compare(-1.1, -1.0) < 0);
  gt_ensure(had_err, gt_double_compare(1+less_than_epsilon, 1.0) == 0);
  gt_ensure(had_err, gt_double_compare(1+less_than_epsilon, -1.0) > 0);
  gt_ensure(had_err, gt_double_compare(-1+less_than_epsilon, -1.0) == 0);
  gt_ensure(had_err, gt_double_compare(-1+less_than_epsilon, 1.0) < 0);

  gt_ensure(had_err, gt_double_smaller_double(1.0, 1.1));
  gt_ensure(had_err, gt_double_smaller_double(-1.0, 1.1));
  gt_ensure(had_err, gt_double_smaller_double(-1.1, -1.0));
  gt_ensure(had_err, !gt_double_smaller_double(-1.0, -1.1));
  gt_ensure(had_err, !gt_double_smaller_double(1.0-less_than_epsilon, 1.0));

  gt_ensure(had_err,  0L == gt_round_to_long(0.5));
  gt_ensure(had_err,  1L == gt_round_to_long(0.51));
  gt_ensure(had_err, -1L == gt_round_to_long(-0.51));
  gt_ensure(had_err,  0L == gt_round_to_long(-0.5));
  gt_ensure(had_err, -2L == gt_round_to_long(-1.5));
  gt_ensure(had_err, -2L == gt_round_to_long(-2.5));
  gt_ensure(had_err, -3L == gt_round_to_long(-2.51));

  return had_err;
}
#endif /* S_SPLINT_S */
