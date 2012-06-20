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

#include "core/ensure.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "extended/elias_gamma.h"

typedef enum {
  LEADING_ZEROS,
  REST
} EliasGammaStatus;

struct GtEliasGammaBitwiseDecoder {
  EliasGammaStatus status;
  unsigned long cur_bit,
                length_in_bits,
                x;
};

GtBittab* gt_elias_gamma_encode(unsigned long x)
{
  GtBittab* code = NULL;
  unsigned long length_in_bits, idx;

  gt_assert(x > 0);
  length_in_bits = (unsigned long) gt_determinebitspervalue(x);
  if (length_in_bits == 0) {
    code = gt_bittab_new(1UL);
    gt_bittab_set_bit(code, 0);
    return code;
  }
  else {
    code = gt_bittab_new(length_in_bits + length_in_bits - 1UL);
    for (idx = 0; idx < length_in_bits; idx++)
      if (((x >> idx) & 1))
        gt_bittab_set_bit(code, idx);
  }
  return code;
}

GtEliasGammaBitwiseDecoder* gt_elias_gamma_bitwise_decoder_new(void)
{
  GtEliasGammaBitwiseDecoder *egbd = gt_malloc(sizeof (*egbd));
  egbd->status = LEADING_ZEROS;
  egbd->cur_bit = 0;
  egbd->length_in_bits = 0;
  egbd->x = 1UL;
  return egbd;
}

void reset_decoder(GtEliasGammaBitwiseDecoder *egbd)
{
  egbd->status = LEADING_ZEROS;
  egbd->cur_bit = 0;
  egbd->x = 1UL;
  egbd->length_in_bits = 0;
}

int gt_elias_gamma_bitwise_decoder_next(GtEliasGammaBitwiseDecoder *egbd,
                                        bool bit, unsigned long *x)
{
  gt_assert(egbd);
  if (egbd->status == LEADING_ZEROS) {
    if (bit == false)
      egbd->length_in_bits += 1;
    else {
      if (egbd->length_in_bits == 0) {
        *x = 1UL;
        reset_decoder(egbd);
        return 0;
      }
      else
        egbd->status = REST;
    }
  }
  else {
    egbd->x = egbd->x << 1;
    if (bit)
      egbd->x = egbd->x | 1;
    egbd->cur_bit++;
    if (egbd->cur_bit == egbd->length_in_bits) {
      *x = egbd->x;
      reset_decoder(egbd);
      return 0;
    }
  }
  return 1;
}

void gt_elias_gamma_bitwise_decoder_delete(GtEliasGammaBitwiseDecoder *egbd)
{
  gt_free(egbd);
}

int gt_elias_gamma_unit_test(GtError *err)
{
  int stat = -1,
      had_err = 0;
  unsigned long idx,
                idx_j,
                unit_test_x_size = 100UL,
                number = unit_test_x_size + 1;
  GtBittab *code;
  GtEliasGammaBitwiseDecoder *egbd = gt_elias_gamma_bitwise_decoder_new();

  for (idx = 1UL; idx <= unit_test_x_size; idx++) {
    code = gt_elias_gamma_encode(idx);
    for (idx_j = gt_bittab_size(code) - 1; idx_j != 0 ; idx_j--) {
      if (gt_bittab_bit_is_set(code, idx_j))
        stat = gt_elias_gamma_bitwise_decoder_next(egbd, true, &number);
      else
        stat = gt_elias_gamma_bitwise_decoder_next(egbd, false, &number);
    }
    if (gt_bittab_bit_is_set(code, idx_j))
      stat = gt_elias_gamma_bitwise_decoder_next(egbd, true, &number);
    else
      stat = gt_elias_gamma_bitwise_decoder_next(egbd, false, &number);

    gt_ensure(had_err, stat == 0);
    gt_ensure(had_err, number == idx);
    gt_bittab_delete(code);
  }
  gt_elias_gamma_bitwise_decoder_delete(egbd);
  return had_err;
}
