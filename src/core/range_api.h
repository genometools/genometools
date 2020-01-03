/*
  Copyright (c) 2005-2009 Gordon Gremme <gordon@gremme.org>
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

#ifndef RANGE_API_H
#define RANGE_API_H

#include <stdbool.h>
#include "core/array_api.h"
#include "core/error_api.h"
#include "core/file_api.h"
#include "core/range_api.h"

/* The <GtRange> class is used to represent genomic ranges in __GenomeTools__.
   Thereby, the <start> must ___always___ be smaller or equal than the <end>. */
typedef struct GtRange GtRange;

struct GtRange {
  GtUword start,
          end;
};

/* Compare <range_a> and <range_b>. Returns 0 if <range_a> equals <range_b>, -1
   if <range_a> starts before <range_b> or (for equal starts) <range_a> ends
   before <range_b>, and 1 else. */
int     gt_range_compare(const GtRange *range_a, const GtRange *range_b);
/* Compare <range_a> and <range_b> with given <delta>.
   Returns 0 if <range_a> equals <range_b> modulo <delta> (i.e., the start and
   end points of <range_a> and <range_b> are at most <delta> bases apart), -1
   if <range_a> starts before <range_b> or (for equal starts) <range_a> ends
   before <range_b>, and 1 else. */
int     gt_range_compare_with_delta(const GtRange *range_a,
                                    const GtRange *range_b,
                                    GtUword delta);
/* Returns <true> if <range_a> and <range_b> overlap, <false> otherwise. */
bool    gt_range_overlap(const GtRange *range_a, const GtRange *range_b);
/* Returns <true> if <range_a> and <range_b> overlap ___at least___ <delta> many
   positions, <false> otherwise. */
bool    gt_range_overlap_delta(const GtRange *range_a,
                               const GtRange *range_b,
                               GtUword delta);
/* Returns <true> if <range_b> is contained in <range_a>, <false> otherwise. */
bool    gt_range_contains(const GtRange *range_a, const GtRange *range_b);
/* Returns <true> if <point> lies within <range>, <false> otherwise. */
bool    gt_range_within(const GtRange *range, GtUword point);
/* Join <range_a> and <range_b> and return the result. */
GtRange gt_range_join(const GtRange *range_a, const GtRange *range_b);
/* Transform start and end of <range> by <offset> and return the result. */
GtRange gt_range_offset(const GtRange *range, GtWord offset);
/* Returns the length of the given <range>. */
GtUword gt_range_length(const GtRange *range);
/* Reorder <range>. */
GtRange gt_range_reorder(GtRange range);
/* Sort an array <ranges> of ranges. */
void    gt_ranges_sort(GtArray *ranges);
/* Sort an array <ranges> of ranges by length. */
void    gt_ranges_sort_by_length_stable(GtArray *ranges);
/* RetuUE if the ranges in <ranges> are sorted. */
bool    gt_ranges_are_sorted(const GtArray *ranges);
/* Returns TRUE if the ranges in <ranges> do not overlap. */
bool    gt_ranges_do_not_overlap(const GtArray *ranges);
/* Retuns TRUE if the ranges in <ranges> are sorted and do not overlap. */
bool    gt_ranges_are_sorted_and_do_not_overlap(const GtArray *ranges);
/* Returns TRUE if the ranges in <ranges_a> and <ranges_b> are equal. */
bool    gt_ranges_are_equal(const GtArray *ranges_a,
                                  const GtArray *ranges_b);
/* Takes a sorted array of ranges and runs the equivalent of uniq on it. */
void    gt_ranges_uniq(GtArray*, const GtArray*);
/* Similar to the previous function, just in place. */
void    gt_ranges_uniq_in_place(GtArray*);
/* Similar to gt_ranges_uniq(), additionally returns an array which contains the
   counts of the occurrences of each elem in the original array. */
GtArray* gt_ranges_uniq_count(GtArray*, const GtArray*);
/* Similar to the previous function, just in place.S */
GtArray* gt_ranges_uniq_in_place_count(GtArray*);
/* Returns TRUE if the ranges in <ranges> are consecutive. */
bool    gt_ranges_are_consecutive(const GtArray *ranges);
/* Returns the sum of the length of the ranges in <ranges>. */
GtUword gt_ranges_total_length(const GtArray *ranges);
/* Returns the length of the boundaries of the ranges in <ranges>. */
GtUword gt_ranges_spanned_length(const GtArray *ranges);
/* Copies ranges <inranges> to the 'opposite' strand, the result being in
   <outranges>. */
void    gt_ranges_copy_to_opposite_strand(GtArray *outranges,
                                          const GtArray *inranges,
                                          GtUword gen_total_length,
                                          GtUword gen_offset);
/* Returns TRUE if all ranges in <ranges> are within the <region>. */
bool    gt_ranges_borders_are_in_region(GtArray *ranges,
                                        const GtRange *region);
/* Prints a representation of the ranges <ranges> to <outfp>. */
void    gt_ranges_show(GtArray *ranges, GtFile *outfp);

int     gt_range_unit_test(GtError*);

#endif
