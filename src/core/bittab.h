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

typedef struct GT_Bittab GT_Bittab;

GT_Bittab*    gt_bittab_new(unsigned long num_of_bits);
void          gt_bittab_set_bit(GT_Bittab*, unsigned long);
void          gt_bittab_unset_bit(GT_Bittab*, unsigned long);
/* a = ~b */
void          gt_bittab_complement(GT_Bittab*, const GT_Bittab*);
/* a = b */
void          gt_bittab_equal(GT_Bittab*, const GT_Bittab*);
/* a = b & c */
void          gt_bittab_and(GT_Bittab*, const GT_Bittab*, const GT_Bittab*);
/* a = b | c */
void          gt_bittab_or(GT_Bittab*, const GT_Bittab*, const GT_Bittab*);
/* a = b & ~c */
void          gt_bittab_nand(GT_Bittab*, const GT_Bittab*, const GT_Bittab*);
/* a &= b */
void          gt_bittab_and_equal(GT_Bittab*, const GT_Bittab*);
/* a |= b */
void          gt_bittab_or_equal(GT_Bittab*, const GT_Bittab*);
/* a <<= 1 */
void          gt_bittab_shift_left_equal(GT_Bittab*);
/* a >>= 1 */
void          gt_bittab_shift_right_equal(GT_Bittab*);
/* a = 0 */
void          gt_bittab_unset(GT_Bittab*);
void          gt_bittab_show(const GT_Bittab*, FILE*);
void          gt_bittab_get_all_bitnums(const GT_Bittab*, GtArray*);
bool          gt_bittab_bit_is_set(const GT_Bittab*, unsigned long);
bool          gt_bittab_is_true(const GT_Bittab*);
bool          gt_bittab_cmp(const GT_Bittab*, const GT_Bittab*);
unsigned long gt_bittab_get_first_bitnum(const GT_Bittab*);
unsigned long gt_bittab_get_last_bitnum(const GT_Bittab*);
unsigned long gt_bittab_get_next_bitnum(const GT_Bittab*, unsigned long);
unsigned long gt_bittab_count_set_bits(GT_Bittab*);
unsigned long gt_bittab_size(GT_Bittab*);
int           gt_bittab_example(GT_Error*);
int           gt_bittab_unit_test(GT_Error*);
void          gt_bittab_delete(GT_Bittab*);

#endif
