/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef BITTAB_H
#define BITTAB_H

#include <stdbool.h>
#include <stdio.h>
#include "core/array.h"

typedef struct GtBittab GtBittab;

GtBittab*    gt_bittab_new(unsigned long num_of_bits);
void          gt_bittab_set_bit(GtBittab*, unsigned long);
void          gt_bittab_unset_bit(GtBittab*, unsigned long);
/* a = ~b */
void          gt_bittab_complement(GtBittab*, const GtBittab*);
/* a = b */
void          gt_bittab_equal(GtBittab*, const GtBittab*);
/* a = b & c */
void          gt_bittab_and(GtBittab*, const GtBittab*, const GtBittab*);
/* a = b | c */
void          gt_bittab_or(GtBittab*, const GtBittab*, const GtBittab*);
/* a = b & ~c */
void          gt_bittab_nand(GtBittab*, const GtBittab*, const GtBittab*);
/* a &= b */
void          gt_bittab_and_equal(GtBittab*, const GtBittab*);
/* a |= b */
void          gt_bittab_or_equal(GtBittab*, const GtBittab*);
/* a <<= 1 */
void          gt_bittab_shift_left_equal(GtBittab*);
/* a >>= 1 */
void          gt_bittab_shift_right_equal(GtBittab*);
/* a = 0 */
void          gt_bittab_unset(GtBittab*);
void          gt_bittab_show(const GtBittab*, FILE*);
void          gt_bittab_get_all_bitnums(const GtBittab*, GtArray*);
bool          gt_bittab_bit_is_set(const GtBittab*, unsigned long);
bool          gt_bittab_is_true(const GtBittab*);
bool          gt_bittab_cmp(const GtBittab*, const GtBittab*);
unsigned long gt_bittab_get_first_bitnum(const GtBittab*);
unsigned long gt_bittab_get_last_bitnum(const GtBittab*);
unsigned long gt_bittab_get_next_bitnum(const GtBittab*, unsigned long);
unsigned long gt_bittab_count_set_bits(const GtBittab*);
unsigned long gt_bittab_size(GtBittab*);
int           gt_bittab_example(GtError*);
int           gt_bittab_unit_test(GtError*);
void          gt_bittab_delete(GtBittab*);

#endif
