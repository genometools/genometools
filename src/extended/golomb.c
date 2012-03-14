/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>

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
#include "core/divmodmul.h"
#include "core/ensure.h"
#include "core/ma_api.h"
#include "extended/golomb.h"

struct GtGolomb {
  unsigned long m;
};

typedef enum {
  IN_Q,
  REMAINDER
} GolombStatus;

struct GtGolombBitwiseDecoder {
  unsigned long m,
                cur_r_bit,
                q,
                r,
                l;
  GolombStatus status;
};

GtGolomb* gt_golomb_new(unsigned long m)
{
  GtGolomb *gc;
  gt_assert(m > 0);
  gc = gt_malloc(sizeof (*gc));
  gc->m = m;
  return gc;
}

unsigned long gt_golomb_get_m(GtGolomb *gc)
{
  gt_assert(gc);
  return gc->m;
}

GtBittab* gt_golomb_encode(GtGolomb *gc, unsigned long x)
{
  unsigned long q,
                r,
                l,
                bit;

  int idx_i;
  GtBittab *code;
  gt_assert(gc);

  q = (unsigned long) floor (x / gc->m);
  r = x - q * gc->m;

  l = (unsigned long) ceil(log(gc->m) * 1/log(2));

  if (l == 0) {
    code = gt_bittab_new((q + 1) + 1);
  }
  else if (r < (GT_POW2(l)) - gc->m) {
    code = gt_bittab_new((q + 1) + (l - 1));

    for (idx_i = 0; idx_i < l - 1; idx_i++) {
      bit = ((r >> idx_i) & 1);
      if (bit == 1)
        gt_bittab_set_bit(code, idx_i);
    }
  }
  else {
    code = gt_bittab_new((q + 1) + l);
    r = r + (GT_POW2(l)) - gc->m;
    for (idx_i = 0; idx_i < l; idx_i++) {
      bit = ((r >> idx_i) & 1);
      if (bit == 1)
        gt_bittab_set_bit(code, idx_i);
    }
  }

  for (idx_i = gt_bittab_size(code) - 1;
       idx_i >= gt_bittab_size(code) - q;
       idx_i--)
    gt_bittab_set_bit(code, idx_i);
  return code;
}

GtGolombBitwiseDecoder *gt_golomb_bitwise_decoder_new(GtGolomb *gc)
{
  GtGolombBitwiseDecoder *gbwd;
  gt_assert(gc);
  gbwd = gt_malloc(sizeof (*gbwd));
  gbwd->m = gc->m;
  gbwd->status = IN_Q;
  gbwd->cur_r_bit = 0;
  gbwd->q = 0;
  gbwd->r = 0;
  gbwd->l = (unsigned long) ceil(log(gc->m) * 1/log(2));
  return gbwd;
}

int gt_golomb_bitwise_decoder_next(GtGolombBitwiseDecoder *gbwd,
                                   bool bit, unsigned long *x)
{
  gt_assert(gbwd);
  if (gbwd->status == IN_Q) {
    if (bit == true)
      gbwd->q += 1;
    else
      gbwd->status = REMAINDER;
    return 1;
  }
  else if (gbwd->status == REMAINDER) {
    gbwd->r = gbwd->r << 1;
    if (bit == true)
      gbwd->r = gbwd->r | 1;
    gbwd->cur_r_bit++;

    if (gbwd->l == 0) {
      *x = gbwd->q * gbwd->m + gbwd->r;
      gbwd->status = IN_Q;
      gbwd->q = 0;
      gbwd->r = 0;
      gbwd->cur_r_bit = 0;
      return 0;
    }
    else if (gbwd->cur_r_bit == gbwd->l - 1) {
      if (gbwd->r < GT_POW2(gbwd->l) - gbwd->m) {
        *x = gbwd->q * gbwd->m + gbwd->r;
        gbwd->status = IN_Q;
        gbwd->q = 0;
        gbwd->r = 0;
        gbwd->cur_r_bit = 0;
        return 0;
      }
      else
        return 1;
    }
    else if (gbwd->cur_r_bit == gbwd->l) {
      gbwd->r -= ((GT_POW2(gbwd->l)) - gbwd->m);
      *x = gbwd->q * gbwd->m + gbwd->r;
      gbwd->status = IN_Q;
      gbwd->q = 0;
      gbwd->r = 0;
      gbwd->cur_r_bit = 0;
      return 0;
    }
    return 1;
  }
  return 1;
}

void gt_golomb_delete(GtGolomb *gc)
{
  if (!gc)
    return;
  gt_free(gc);
}

void gt_golomb_bitwise_decoder_delete(GtGolombBitwiseDecoder *gbwd)
{
  if (!gbwd)
    return;
  gt_free(gbwd);
}

int gt_golomb_unit_test(GtError *err)
{
  int had_err = 0,
      stat = -1,
      idx_i;
  GtGolomb *gc;
  GtBittab *code = NULL;
  GtGolombBitwiseDecoder *gbwd = NULL;
  unsigned long number,
                idx_j,
                idx_k,
                unit_test_x_size = 100,
                unit_test_b_size = 100;

  for (idx_j = 1; idx_j <= unit_test_b_size; idx_j++) {
    gc = gt_golomb_new(idx_j);
    gbwd = gt_golomb_bitwise_decoder_new(gc);
    for (idx_k = 0; idx_k <= unit_test_x_size; idx_k++) {
      code = gt_golomb_encode(gc, idx_k);
      for (idx_i = gt_bittab_size(code) - 1; idx_i >=  0 ; idx_i--) {
        if (gt_bittab_bit_is_set(code, idx_i))
          stat = gt_golomb_bitwise_decoder_next(gbwd, true, &number);
        else
          stat = gt_golomb_bitwise_decoder_next(gbwd, false, &number);
      }
      gt_ensure(had_err, stat == 0);
      gt_ensure(had_err, number == idx_k);
      gt_bittab_delete(code);
    }
    gt_golomb_bitwise_decoder_delete(gbwd);
    gt_golomb_delete(gc);
  }
  return had_err;
}
