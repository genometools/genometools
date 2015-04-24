/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef MATHSUPPORT_H
#define MATHSUPPORT_H

#include <stdbool.h>
#include <math.h>
#include "core/error_api.h"
#include "core/types_api.h"

/* Returns the log of the sum of two log probabilities. */
double       gt_logsum(double p1, double p2);
bool         gt_double_equals_one(double);
bool         gt_double_equals_double(double, double);
int          gt_double_compare(double, double);
bool         gt_double_smaller_double(double, double);
bool         gt_double_larger_double(double, double);

/* Returns a random number between 0 and maximal_value. */
GtUword      gt_rand_max(GtUword maximal_value);
/* Returns a random double between 0.0 and maximal_value. */
double       gt_rand_max_double(double maximal_value);
/* Returns a random double between 0.0 and 1.0. */
double       gt_rand_0_to_1(void);
/* Returns a random character from 'a' to 'z'. */
char         gt_rand_char(void);
/* Retuns the log base 2 of an integer <maxvalue> in O(wordsize) operations */
unsigned int gt_determinebitspervalue(GtUword maxvalue);
/* Determine pow(base,exponent) for small values of exponent */
GtUword      gt_power_for_small_exponents(unsigned int base,
                                          unsigned int exponent);
/* Return <x> rounded to the nearest long integer, similar to
   <lrint()> which may not be available on older glibc versions. */
GtWord       gt_round_to_long(double x);
/* Compute the greatest common divisor of two unsigned integers */
unsigned int gt_gcd_uint(unsigned int m, unsigned int n);
/* Compute the least common multiplier of two unsigned integers */
unsigned int gt_lcm_uint(unsigned int m, unsigned int n);
/* Compute the logarith of <x> to the base <b> */
double       gt_log_base(double x, double b);

int          gt_mathsupport_unit_test(GtError *err);

#endif
