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
                l,
                x;
};

GtBittab* gt_elias_gamma_encode(unsigned long x)
{
  GtBittab* code = NULL;
  unsigned int l;
  long i;

  gt_assert(x > 0);
  l = gt_determinebitspervalue(x);
  if (l == 0) {
    code = gt_bittab_new(1);
    gt_bittab_set_bit(code, 0);
    return code;
  }
  else {
    code = gt_bittab_new(l + l - 1);
    for (i = 0; i < l; i++)
      if (((x >> i) & 1))
        gt_bittab_set_bit(code, i);
  }
  return code;
}

GtEliasGammaBitwiseDecoder* gt_elias_gamma_bitwise_decoder_new()
{
  GtEliasGammaBitwiseDecoder *egbd = gt_malloc(sizeof (*egbd));
  egbd->status = LEADING_ZEROS;
  egbd->cur_bit = 0;
  egbd->l = 0;
  egbd->x = 1;
  return egbd;
}

int gt_elias_gamma_bitwise_decoder_next(GtEliasGammaBitwiseDecoder *egbd,
                                        bool bit, unsigned long *x)
{
  gt_assert(egbd);
  if (egbd->status == LEADING_ZEROS) {
    if (bit == false)
      egbd->l += 1;
    else {
      if (egbd->l == 0) {
        *x = 1;
        egbd->status = LEADING_ZEROS;
        egbd->cur_bit = 0;
        egbd->x = 1;
        egbd->l = 0;
        return 0;
      }
      else
        egbd->status = REST;
    }
    return 1;
  }
  else {
    egbd->x = egbd->x << 1;
    if (bit == true)
      egbd->x = egbd->x | 1;
    egbd->cur_bit++;
    if (egbd->cur_bit == egbd->l) {
      *x = egbd->x;
      egbd->status = LEADING_ZEROS;
      egbd->cur_bit = 0;
      egbd->x = 1;
      egbd->l = 0;
      return 0;
    }
    else
      return 1;
  }
  return 0;
}

void gt_elias_gamma_bitwise_decoder_delete(GtEliasGammaBitwiseDecoder *egbd)
{
  if (!egbd)
    return;
  gt_free(egbd);
}

int gt_elias_gamma_unit_test(GtError *err)
{
  long i,
       j;
  unsigned long number,
                unit_test_x_size = 100;
  int stat = -1,
      had_err = 0;
  GtBittab *code;
  GtEliasGammaBitwiseDecoder *egbd = gt_elias_gamma_bitwise_decoder_new();
  for (i = 1; i <= unit_test_x_size; i++) {
    code = gt_elias_gamma_encode(i);
    for (j = gt_bittab_size(code) - 1; j >=  0 ; j--) {
      if (gt_bittab_bit_is_set(code, j))
        stat = gt_elias_gamma_bitwise_decoder_next(egbd, true, &number);
      else
        stat =
              gt_elias_gamma_bitwise_decoder_next(egbd, false, &number);
    }
    gt_ensure(had_err, stat == 0);
    gt_ensure(had_err, number == i);
    gt_bittab_delete(code);
  }
  gt_elias_gamma_bitwise_decoder_delete(egbd);
  return had_err;
}
