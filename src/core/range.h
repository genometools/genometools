/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef RANGE_H
#define RANGE_H

#include <stdbool.h>
#include "core/array.h"

typedef struct {
  unsigned long start,
                end;
} GT_Range;

int           range_compare(GT_Range, GT_Range);
int           range_compare_ptr(const GT_Range*, const GT_Range*);
int           range_compare_with_delta(GT_Range, GT_Range, unsigned long delta);
bool          range_overlap(GT_Range, GT_Range);
bool          range_contains(GT_Range, GT_Range);
bool          range_within(GT_Range, unsigned long);
GT_Range         range_join(GT_Range, GT_Range);
GT_Range         range_offset(GT_Range, long offset);
GT_Range         range_reorder(GT_Range);
unsigned long range_length(GT_Range);
int           range_unit_test(GT_Error*);

void          ranges_sort(GT_Array*);
void          ranges_sort_by_length_stable(GT_Array*);
bool          ranges_are_sorted(const GT_Array*);
bool          ranges_do_not_overlap(const GT_Array*);
bool          ranges_are_sorted_and_do_not_overlap(const GT_Array*);
bool          ranges_are_equal(const GT_Array*, const GT_Array*);

/* takes a sorted array of ranges and runs the equivalent of uniq on it. The
   result is returned. */
void          ranges_uniq(GT_Array*, const GT_Array*);
/* similar to the previous function, just in place */
void          ranges_uniq_in_place(GT_Array*);
/* similar to ranges_uniq(), additionally returns an array which contains the
   counts of the occurrences of each elem in the original array */
GT_Array*     ranges_uniq_count(GT_Array*, const GT_Array*);
/* similar to the previous function, just in place */
GT_Array*     ranges_uniq_in_place_count(GT_Array*);

#endif
