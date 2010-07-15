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

#ifndef BITTAB_API_H
#define BITTAB_API_H

#include <stdbool.h>
#include <stdio.h>
#include "core/array_api.h"

/* Implements arbitrary-length bit arrays and various operations on them. */
typedef struct GtBittab GtBittab;

/* Creates a new <GtBittab> of length <num_of_bits>, initialised to 0 */
GtBittab*     gt_bittab_new(unsigned long num_of_bits);

/* Sets bit <i> in <b> to 1. */
void          gt_bittab_set_bit(GtBittab *b, unsigned long i);

/* Sets bit <i> in <b> to 0. */
void          gt_bittab_unset_bit(GtBittab *b, unsigned long i);

/* Sets <a> to be the complement of <b>. */
void          gt_bittab_complement(GtBittab *a, const GtBittab *b);

/* Sets <a> to be equal to <b>. */
void          gt_bittab_equal(GtBittab *a, const GtBittab *b);

/* Sets <a> to be the bitwise AND of <b> and <c>. */
void          gt_bittab_and(GtBittab *a, const GtBittab *b, const GtBittab *c);

/* Sets <a> to be the bitwise OR of <b> and <c>. */
void          gt_bittab_or(GtBittab *a, const GtBittab *b, const GtBittab *c);

/* Sets <a> to be <b> NAND <c>. */
void          gt_bittab_nand(GtBittab *a, const GtBittab *b, const GtBittab *c);

/* Sets <a> to be the bitwise AND of <a> and <b>. */
void          gt_bittab_and_equal(GtBittab *a, const GtBittab *b);

/* Sets <a> to be the bitwise OR of <a> and <b>. */
void          gt_bittab_or_equal(GtBittab *a, const GtBittab *b);

/* Shifts <b> by one position to the left. */
void          gt_bittab_shift_left_equal(GtBittab *b);

/* Shifts <b> by one position to the right. */
void          gt_bittab_shift_right_equal(GtBittab *b);

/* Sets all bits in <b> to 0. */
void          gt_bittab_unset(GtBittab *b);

/* Outputs a representation of <b> to <fp>. */
void          gt_bittab_show(const GtBittab *b, FILE *fp);

/* Fills <a> with the indices of all set bits in <b>. */
void          gt_bittab_get_all_bitnums(const GtBittab *b, GtArray *a);

/* Returns <true> if bit <i> is set in <b>. */
bool          gt_bittab_bit_is_set(const GtBittab *b, unsigned long i);

/* Returns <true> if bittabs <a> and <b> are identical. */
bool          gt_bittab_cmp(const GtBittab *a, const GtBittab *b);

/* Returns the index of the first set bit in <b>. */
unsigned long gt_bittab_get_first_bitnum(const GtBittab *b);

/* Returns the index of the last set bit in <b>. */
unsigned long gt_bittab_get_last_bitnum(const GtBittab *b);

/* Returns the index of the next set bit in <b> with an index greater
   than <i>. */
unsigned long gt_bittab_get_next_bitnum(const GtBittab *b, unsigned long i);

/* Returns the number of set bits in <b>. */
unsigned long gt_bittab_count_set_bits(const GtBittab *b);

/* Returns the total number of bits of <b>. */
unsigned long gt_bittab_size(GtBittab *b);

/* Deletes <b> and frees all allocated space. */
void          gt_bittab_delete(GtBittab *b);

#endif
