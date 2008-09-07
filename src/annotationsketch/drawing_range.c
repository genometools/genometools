/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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

#include <assert.h>
#include "core/mathsupport.h"
#include "annotationsketch/drawing_range.h"

int drawing_range_compare(DrawingRange gt_range_a, DrawingRange gt_range_b)
{
  assert(gt_range_a.start <= gt_range_a.end && gt_range_b.start <= gt_range_b.end);

  if (double_equals_double(gt_range_a.start, gt_range_b.start)
        && double_equals_double(gt_range_a.end, gt_range_b.end))
    return 0; /* gt_range_a == gt_range_b */

  if ((gt_range_a.start < gt_range_b.start) ||
      (double_equals_double(gt_range_a.start, gt_range_b.start)
         && (gt_range_a.end < gt_range_b.end)))
    return -1; /* gt_range_a < gt_range_b */

  return 1; /* gt_range_a > gt_range_b */
}

bool drawing_range_overlap(DrawingRange gt_range_a, DrawingRange gt_range_b)
{
  if (gt_range_a.start <= gt_range_b.end && gt_range_a.end >= gt_range_b.start)
    return true;
  return false;
}

bool drawing_range_contains(DrawingRange gt_range_a, DrawingRange gt_range_b)
{
  if (gt_range_a.start <= gt_range_b.start && gt_range_a.end >= gt_range_b.end)
    return true;
  return false;
}

double drawing_range_length(DrawingRange range)
{
  return range.end - range.start;
}
