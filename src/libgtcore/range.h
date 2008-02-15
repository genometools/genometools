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
#include "libgtcore/array.h"

typedef struct {
  unsigned long start,
                end;
} Range;

int           range_compare(Range, Range);
int           range_compare_ptr(const Range*, const Range*);
int           range_compare_with_delta(Range, Range, unsigned long delta);
bool          range_overlap(Range, Range);
bool          range_contains(Range, Range);
bool          range_within(Range, unsigned long);
Range         range_join(Range, Range);
Range         range_offset(Range, long offset);
Range         range_reorder(Range);
unsigned long range_length(Range);
int           range_unit_test(Error*);

void          ranges_sort(Array*);
void          ranges_sort_by_length_stable(Array*);
bool          ranges_are_sorted(const Array*);
bool          ranges_do_not_overlap(const Array*);
bool          ranges_are_sorted_and_do_not_overlap(const Array*);
bool          ranges_are_equal(const Array*, const Array*);

/* takes a sorted array of ranges and runs the equivalent of uniq on it. The
   result is returned. */
void          ranges_uniq(Array*, const Array*);
/* similar to the previous function, just in place */
void          ranges_uniq_in_place(Array*);
/* similar to ranges_uniq(), additionally returns an array which contains the
   counts of the occurrences of each elem in the original array */
Array*        ranges_uniq_count(Array*, const Array*);
/* similar to the previous function, just in place */
Array*        ranges_uniq_in_place_count(Array*);

#endif
