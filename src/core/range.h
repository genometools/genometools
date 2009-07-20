/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/error_api.h"
#include "core/file.h"
#include "core/range_api.h"

GtRange       gt_range_reorder(GtRange);
int           gt_range_unit_test(GtError*);

void          gt_ranges_sort(GtArray*);
void          gt_ranges_sort_by_length_stable(GtArray*);
bool          gt_ranges_are_sorted(const GtArray*);
bool          gt_ranges_do_not_overlap(const GtArray*);
bool          gt_ranges_are_sorted_and_do_not_overlap(const GtArray*);
bool          gt_ranges_are_equal(const GtArray*, const GtArray*);

/* takes a sorted array of ranges and runs the equivalent of uniq on it. The
   result is returned. */
void          gt_ranges_uniq(GtArray*, const GtArray*);
/* similar to the previous function, just in place */
void          gt_ranges_uniq_in_place(GtArray*);
/* similar to gt_ranges_uniq(), additionally returns an array which contains the
   counts of the occurrences of each elem in the original array */
GtArray*      gt_ranges_uniq_count(GtArray*, const GtArray*);
/* similar to the previous function, just in place */
GtArray*      gt_ranges_uniq_in_place_count(GtArray*);

bool          gt_ranges_are_consecutive(const GtArray*);
unsigned long gt_ranges_total_length(const GtArray*);
unsigned long gt_ranges_spanned_length(const GtArray*);
void          gt_ranges_copy_to_opposite_strand(GtArray *outranges,
                                                const GtArray *inranges,
                                                unsigned long gen_total_length,
                                                unsigned long gen_offset);
bool          gt_ranges_borders_are_in_region(GtArray *ranges,
                                              const GtRange *region);
void          gt_ranges_show(GtArray *ranges, GtFile *outfp);

#endif
