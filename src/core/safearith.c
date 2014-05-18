/*
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/ensure.h"
#include "core/safearith.h"
#include "core/types_api.h"
#include "core/unused_api.h"

void gt_safe_default_overflow_handler(const char *src_file, int src_line,
                                      GT_UNUSED void *data)
{
  fprintf(stderr, "%s, l.%d: overflow in operation\n", src_file, src_line);
  exit(EXIT_FAILURE);
}

int gt_safe_abs_check_func(int j, const char *src_file, int src_line,
                      GtOverflowHandlerFunc handler_func, void *data)
{
  int rval = j < 0 ? -j : j;
  if (rval < 0) { /* overflow */
    handler_func(src_file, src_line, data);
  }
  return rval;
}

GtWord gt_safe_labs_check_func(GtWord j, const char *src_file, int src_line,
                             GtOverflowHandlerFunc handler_func, void *data)
{
  GtWord rval = j < 0 ? -j : j;
  if (rval < 0) { /* overflow */
    handler_func(src_file, src_line, data);
  }
  return rval;
}

GtInt64 gt_safe_llabs_check_func(GtInt64 j, const char *src_file, int src_line,
                                 GtOverflowHandlerFunc handler_func, void *data)
{
  GtInt64 rval = j < 0 ? -j : j;
  if (rval < 0) { /* overflow */
    handler_func(src_file, src_line, data);
  }
  return rval;
}

uint32_t gt_safe_mult_u32_check_func(uint32_t a, uint32_t b,
                                     const char *src_file, int src_line,
                                     GtOverflowHandlerFunc handler_func,
                                     void *data)
{
  GtUint64 x = (GtUint64) a * b;
  if (x > 0xffffffff) { /* overflow */
    handler_func(src_file, src_line, data);
  }
  return x;
}

uint64_t gt_safe_mult_u64_check_func(uint64_t a, uint64_t b,
                                     const char *src_file, int src_line,
                                     GtOverflowHandlerFunc handler_func,
                                     void *data)
{
  uint32_t a_hi = a >> 32,
           a_lo = a & 0xffffffff,
           b_hi = b >> 32,
           b_lo = b & 0xffffffff;
  /*
     x = 0xffffffff
     a = a_hi*x + a_lo
     b = b_hi*x + b_lo
     a * b = (a_hi*x + a_lo) * (b_hi*x + b_lo)
           = a_hi*x*b_hi*x + a_hi*x*b_lo + a_lo*b_hi*x + a_lo*b_lo
  */
  if (a_hi && b_hi) {    /* overflow */
    handler_func(src_file, src_line, data);
  }
  a = (uint64_t)(a_hi) * b_lo + (uint64_t)(a_lo) * b_hi;
  if (a > 0xffffffff) {  /* overflow */
    handler_func(src_file, src_line, data);
  }
  return (a << 32) + (uint64_t)(a_lo) * b_lo;
}

GtUword gt_safe_mult_ulong_check_func(GtUword a, GtUword b,
                                            const char *src_file, int src_line,
                                            GtOverflowHandlerFunc handler_func,
                                            void *data)
{
  gt_assert(sizeof (GtUword) == 4 || sizeof (GtUword) == 8);
  if (sizeof (GtUword) == 4) {
    return gt_safe_mult_u32_check_func(a, b, src_file, src_line, handler_func,
                                       data);
  } else { /* sizeof (GtUword) == 8 */
    return gt_safe_mult_u64_check_func(a, b, src_file, src_line, handler_func,
                                       data);
  }
}

GtWord gt_safe_cast2long_check_func(GtUword value, const char *src_file,
                                  int src_line,
                                  GtOverflowHandlerFunc handler_func,
                                  void *data)
{
  if (value > (~0UL >> 1)) {
    handler_func(src_file, src_line, data);
  }
  return value;
}

GtUword gt_safe_cast2ulong_check_func(GtWord value, const char *src_file,
                                            int src_line,
                                            GtOverflowHandlerFunc handler_func,
                                            void *data)
{
  if (value < 0) {
    handler_func(src_file, src_line, data);
  }
  return value;
}

GtUword gt_safe_cast2ulong_64_check_func(uint64_t value,
                                             const char *src_file,
                                             int src_line,
                                             GtOverflowHandlerFunc handler_func,
                                             void *data)
{
  if (value > (uint64_t) ULONG_MAX) {
    handler_func(src_file, src_line, data);
  }
  return (GtUword) value;
}

int gt_safearith_example(GT_UNUSED GtError *err)
{
  GtUword ulong;
  GtWord slong;
  unsigned int a, b, c;
  gt_error_check(err);

  /* safe assignments */
  slong = 256;
  gt_safe_assign(ulong, slong);
  gt_assert(ulong = 256);

  /* safe additions */
  a = 256;
  b = 1;
  gt_safe_add(c, a, b);
  gt_assert(c == 257);

  /* safe subtractions */
  a = 256;
  b = 1;
  gt_safe_sub(c, a, b);
  gt_assert(c == 255);

  /* safe absolutes */
  gt_assert(gt_safe_abs(-256) == 256);

  return 0;
}

int gt_safearith_unit_test(GtError *err)
{
  int had_err = 0;
  gt_error_check(err);

  {
    /* This is not always true, e.g. on PPC and ARM!
       cf. http://www.network-theory.co.uk/docs/gccintro/gccintro_71.html
    gt_ensure(__MIN(char) == -128);
    gt_ensure(__MAX(char) == 127); */
    gt_ensure(__MIN(unsigned char) == 0);
    gt_ensure(__MAX(unsigned char) == 255);

    gt_ensure(__MIN(short) == SHRT_MIN);
    gt_ensure(__MAX(short) == SHRT_MAX);
    gt_ensure(__MIN(unsigned short) == 0);
    gt_ensure(__MAX(unsigned short) == USHRT_MAX);

    gt_ensure(__MIN(int) == INT_MIN);
    gt_ensure(__MAX(int) == INT_MAX);
    gt_ensure(__MIN(unsigned int) == 0);
    gt_ensure(__MAX(unsigned int) == UINT_MAX);

    gt_ensure(__MIN(GtWord) == LONG_MIN);
    gt_ensure(__MAX(GtWord) == LONG_MAX);
    gt_ensure(__MIN(GtUword) == 0);
    gt_ensure(__MAX(GtUword) == ULONG_MAX);

#ifdef LLONG_MIN
    gt_ensure(__MIN(GtInt64) == LLONG_MIN);
    gt_ensure(__MAX(GtInt64) == LLONG_MAX);
    gt_ensure(__MIN(GtUint64) == 0);
    gt_ensure(__MAX(GtUint64) == ULLONG_MAX);
#endif
  }

  {
    GtUword ulong;
    GtWord slong;

    slong = -1;
    gt_ensure(assign(ulong, slong));

    ulong = 0;
    slong = LONG_MAX;
    gt_ensure(!assign(ulong, slong) && ulong == LONG_MAX);

    ulong = ULONG_MAX;
    gt_ensure(assign(slong, ulong));

    slong = 0;
    ulong = LONG_MAX;
    gt_ensure(!assign(slong, ulong) && slong == LONG_MAX);
  }

  {
    int x;
    gt_ensure(add_of(x, INT_MAX, 1));
    gt_ensure(add_of(x, INT_MAX, 256));
    gt_ensure(add_of(x, INT_MAX, INT_MAX));

    x = 0; gt_ensure(!add_of(x, INT_MAX - 1, 1) && x == INT_MAX);
    x = 0; gt_ensure(!add_of(x, INT_MAX - 256, 256) && x == INT_MAX);
    x = 0; gt_ensure(!add_of(x, INT_MAX, 0) && x == INT_MAX);

    gt_ensure(add_of(x, 0x100000000ll, 0x100000000ll));
    x = INT_MAX;
    gt_ensure(!add_of(x, 0x100000000ll, -0x100000000ll) && x == 0);

    gt_ensure(sub_of(x, INT_MIN, 1));
    gt_ensure(sub_of(x, INT_MIN, 256));
    gt_ensure(sub_of(x, INT_MIN, INT_MAX));

    x = 0; gt_ensure(!sub_of(x, INT_MIN + 1, 1) && x == INT_MIN);
    x = 0; gt_ensure(!sub_of(x, INT_MIN + 256, 256) && x == INT_MIN);
    x = 0; gt_ensure(!sub_of(x, INT_MIN, 0) && x == INT_MIN);
  }

  {
    unsigned int x;
    gt_ensure(add_of(x, UINT_MAX, 1));
    gt_ensure(add_of(x, UINT_MAX, 256));
    gt_ensure(add_of(x, UINT_MAX, UINT_MAX));

    x = 0; gt_ensure(!add_of(x, UINT_MAX - 1, 1) && x == UINT_MAX);
    x = 0; gt_ensure(!add_of(x, UINT_MAX - 256, 256) && x == UINT_MAX);
    x = 0; gt_ensure(!add_of(x, UINT_MAX, 0) && x == UINT_MAX);

    gt_ensure(add_of(x, 0x100000000ll, 0x100000000ll));
    x = UINT_MAX;
    gt_ensure(!add_of(x, 0x100000000ll, -0x100000000ll) && x == 0);

    gt_ensure(sub_of(x, 0, 1));
    gt_ensure(sub_of(x, 0, 256));
    gt_ensure(sub_of(x, 0, UINT_MAX));

    x = 0; gt_ensure(!sub_of(x, 1, 1) && x == 0);
    x = 0; gt_ensure(!sub_of(x, 256, 256) && x == 0);
    x = 0; gt_ensure(!sub_of(x, 0, 0) && x == 0);
  }

  {
    int i;
    GtWord l;
    GtInt64 ll;

    i = gt_safe_abs(0);
    gt_ensure(i == 0);

    i = gt_safe_abs(-1);
    gt_ensure(i == 1);

    i = gt_safe_abs(INT_MIN + 1);
    gt_ensure(i == INT_MAX);

    l = gt_safe_labs(0);
    gt_ensure(l == 0);

    l = gt_safe_labs(-1);
    gt_ensure(l == 1);

    l = gt_safe_labs(LONG_MIN + 1);
    gt_ensure(l == LONG_MAX);

    ll = gt_safe_llabs(0);
    gt_ensure(ll == 0);

    ll = gt_safe_llabs(-1);
    gt_ensure(ll == 1);

#ifdef LLONG_MIN
    ll = gt_safe_llabs(LLONG_MIN + 1);
    gt_ensure(ll == LLONG_MAX);
#endif
  }

  return had_err;
}
