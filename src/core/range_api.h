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

#ifndef RANGE_API_H
#define RANGE_API_H

#include <stdbool.h>
#include "core/array_api.h"

typedef struct {
  unsigned long start,
                end;
} GT_Range;

int           gt_range_compare(GT_Range, GT_Range);
int           gt_range_compare_ptr(const GT_Range*, const GT_Range*);
int           gt_range_compare_with_delta(GT_Range, GT_Range,
                                          unsigned long delta);
bool          gt_range_overlap(GT_Range, GT_Range);
bool          gt_range_contains(GT_Range, GT_Range);
bool          gt_range_within(GT_Range, unsigned long);
GT_Range      gt_range_join(GT_Range, GT_Range);
GT_Range      gt_range_offset(GT_Range, long offset);
unsigned long gt_range_length(GT_Range);

#endif
