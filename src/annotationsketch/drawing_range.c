/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/mathsupport.h"
#include "annotationsketch/drawing_range.h"

int gt_drawing_range_compare(GtDrawingRange range_a, GtDrawingRange range_b)
{
  gt_assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  if (gt_double_equals_double(range_a.start, range_b.start) &&
      gt_double_equals_double(range_a.end, range_b.end))
    return 0; /* range_a == range_b */

  if (gt_double_smaller_double(range_a.start, range_b.start) ||
      (gt_double_equals_double(range_a.start, range_b.start)
       && gt_double_smaller_double(range_a.end, range_b.end)))
    return -1; /* range_a < range_b */

  return 1; /* range_a > range_b */
}

bool gt_drawing_range_overlap(GtDrawingRange range_a, GtDrawingRange range_b)
{
  if (range_a.start <= range_b.end && range_a.end >= range_b.start)
    return true;
  return false;
}

bool gt_drawing_range_contains(GtDrawingRange range_a, GtDrawingRange range_b)
{
  if (range_a.start <= range_b.start && range_a.end >= range_b.end)
    return true;
  return false;
}

double gt_drawing_range_length(GtDrawingRange range)
{
  return range.end - range.start;
}
