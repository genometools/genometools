/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>

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

/* TODO write stream reader, that returns only, when complete code is read */
#include <math.h>
#include "core/divmodmul.h"
#include "core/ensure.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "extended/golomb.h"

struct GtGolomb {
  unsigned long median;
  unsigned long len;
  unsigned long two_pow_len;
};

typedef enum {
  IN_Q,
  REMAINDER
} GolombStatus;

struct GtGolombBitwiseDecoder {
  unsigned long median,
                cur_r_bit,
                quotient,
                remain,
                len,
                two_pow_len;
  GolombStatus status;
};

GtGolomb* gt_golomb_new(unsigned long median)
{
  GtGolomb *golomb;
  gt_assert(median > 0);
  golomb = gt_malloc(sizeof (*golomb));
  golomb->median = median;
  golomb->len = (unsigned long) ceil(GT_LOG2((double) golomb->median));
  golomb->two_pow_len = GT_POW2(golomb->len);
  return golomb;
}

unsigned long gt_golomb_get_m(GtGolomb *golomb)
{
  gt_assert(golomb);
  return golomb->median;
}

GtBittab* gt_golomb_encode(GtGolomb *golomb, unsigned long x)
{
  unsigned long quotient,
                remain,
                mask,
                idx_i;
  GtBittab *code;
  gt_assert(golomb);

  quotient = x / golomb->median;
  remain = x - quotient * golomb->median;

  if (golomb->len == 0) {
    code = gt_bittab_new((quotient + 1) + 1);
  }
  else if (remain < golomb->two_pow_len - golomb->median) {
    code = gt_bittab_new((quotient + 1) + (golomb->len - 1));

    mask = 1UL << (golomb->len - 2);
    for (idx_i = 0; idx_i < golomb->len - 1; idx_i++) {
      if (remain & mask)
        gt_bittab_set_bit(code, quotient + 1 + idx_i);
      mask >>= 1;
    }
  }
  else {
    code = gt_bittab_new((quotient + 1) + golomb->len);

    remain = remain + (golomb->two_pow_len) - golomb->median;

    mask = 1UL << (golomb->len - 1);
    for (idx_i = 0; idx_i < golomb->len; idx_i++) {
      if (remain & mask)
        gt_bittab_set_bit(code, quotient + 1 + idx_i);
      mask >>= 1;
    }
  }

  /* write quotient 1s */
  /* for (idx_i = gt_bittab_size(code) - 1;
       idx_i >= gt_bittab_size(code) - quotient;
       idx_i--)
    gt_bittab_set_bit(code, idx_i); */
  for (idx_i = 0; idx_i < quotient; idx_i++)
    gt_bittab_set_bit(code, idx_i);
  return code;
}

GtGolombBitwiseDecoder *gt_golomb_bitwise_decoder_new(GtGolomb *golomb)
{
  GtGolombBitwiseDecoder *gbwd;
  gt_assert(golomb);
  gbwd = gt_malloc(sizeof (*gbwd));
  gbwd->median = golomb->median;
  gbwd->status = IN_Q;
  gbwd->cur_r_bit = 0;
  gbwd->quotient = 0;
  gbwd->remain = 0;
  gbwd->len = golomb->len;
  gbwd->two_pow_len = golomb->two_pow_len;
  return gbwd;
}

int gt_golomb_bitwise_decoder_next(GtGolombBitwiseDecoder *gbwd,
                                   bool bit, unsigned long *x)
{
  gt_assert(gbwd);
  if (gbwd->status == IN_Q) {
    if (bit)
      gbwd->quotient += 1;
    else
      gbwd->status = REMAINDER;
    return 1;
  }
  else if (gbwd->status == REMAINDER) {
    gbwd->remain = gbwd->remain << 1;
    if (bit)
      gbwd->remain = gbwd->remain | 1;
    gbwd->cur_r_bit++;

    if (gbwd->len == 0) {
      *x = gbwd->quotient * gbwd->median + gbwd->remain;
      gbwd->status = IN_Q;
      gbwd->quotient = 0;
      gbwd->remain = 0;
      gbwd->cur_r_bit = 0;
      return 0;
    }
    else if (gbwd->cur_r_bit == gbwd->len - 1) {
      if (gbwd->remain < GT_POW2(gbwd->len) - gbwd->median) {
        *x = gbwd->quotient * gbwd->median + gbwd->remain;
        gbwd->status = IN_Q;
        gbwd->quotient = 0;
        gbwd->remain = 0;
        gbwd->cur_r_bit = 0;
        return 0;
      }
      else
        return 1;
    }
    else if (gbwd->cur_r_bit == gbwd->len) {
      gbwd->remain -= ((GT_POW2(gbwd->len)) - gbwd->median);
      *x = gbwd->quotient * gbwd->median + gbwd->remain;
      gbwd->status = IN_Q;
      gbwd->quotient = 0;
      gbwd->remain = 0;
      gbwd->cur_r_bit = 0;
      return 0;
    }
    return 1;
  }
  return -1;
}

void gt_golomb_delete(GtGolomb *golomb)
{
  gt_free(golomb);
}

void gt_golomb_bitwise_decoder_delete(GtGolombBitwiseDecoder *gbwd)
{
  gt_free(gbwd);
}

int gt_golomb_unit_test(GtError *err)
{
  int had_err = 0,
      stat = -1;
  GtGolomb *golomb;
  GtBittab *code = NULL;
  GtGolombBitwiseDecoder *gbwd = NULL;
  unsigned long unit_test_x_size = 128UL,
                idx_i,
                idx_j,
                idx_k,
                number = unit_test_x_size + 1,
                unit_test_b_size = 256UL;

  for (idx_j = 1UL; !had_err && idx_j <= unit_test_b_size; idx_j++) {
    golomb = gt_golomb_new(idx_j);
    gbwd = gt_golomb_bitwise_decoder_new(golomb);
    for (idx_k = 0; !had_err && idx_k <= unit_test_x_size; idx_k++) {
      code = gt_golomb_encode(golomb, idx_k);
      for (idx_i = 0 ; !had_err && idx_i < gt_bittab_size(code) ; idx_i++) {
        if (gt_bittab_bit_is_set(code, idx_i))
          stat = gt_golomb_bitwise_decoder_next(gbwd, true, &number);
        else
          stat = gt_golomb_bitwise_decoder_next(gbwd, false, &number);
        gt_ensure(had_err, stat != -1);
      }

      gt_ensure(had_err, stat == 0);
      gt_ensure(had_err, number == idx_k);
      gt_bittab_delete(code);
    }
    gt_golomb_bitwise_decoder_delete(gbwd);
    gt_golomb_delete(golomb);
  }
  return had_err;
}
