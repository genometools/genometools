/*
  Copyright (c) 2005-2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef RANGE_H
#define RANGE_H

#include "array.h"

typedef struct {
  unsigned long start,
                end;
} Range;

int           range_compare(Range, Range);
int           range_compare_ptr(const Range*, const Range*);
unsigned int  range_overlap(Range, Range);
unsigned int  range_contains(Range, Range);
Range         range_join(Range, Range);
Range         range_offset(Range, long offset, unsigned long line_number);
unsigned long range_length(Range);
int           range_unit_test(void);

void          ranges_sort(Array*);
void          ranges_sort_by_length_stable(Array*);
unsigned int  ranges_are_sorted(const Array*);
unsigned int  ranges_do_not_overlap(const Array*);
unsigned int  ranges_are_sorted_and_do_not_overlap(const Array*);
unsigned int  ranges_are_equal(const Array*, const Array*);

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
