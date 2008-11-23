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

#include <stdio.h>
#include "annotationsketch/cliptype.h"
#include "annotationsketch/coords.h"
#include "core/minmax.h"

double gt_coords_convert_point(GtRange viewrange, long pos)
{
  return ((double) (((long) pos -(long) viewrange.start)))
                  / ((double) gt_range_length(&viewrange));
}

GtDrawingRange gt_coords_calc_generic_range(GtRange node_range,
                                            GtRange viewrange)
{
  GtDrawingRange converted_range;
  converted_range.clip = CLIPPED_NONE;
  node_range.end++;
  /* scale coordinates to target image width */
  if (node_range.start < viewrange.start )
    converted_range.clip = CLIPPED_LEFT;
  converted_range.start = gt_coords_convert_point(viewrange, node_range.start);

  if (node_range.end > viewrange.end+1)
    converted_range.clip = (converted_range.clip == CLIPPED_LEFT ?
                                                      CLIPPED_BOTH :
                                                      CLIPPED_RIGHT);
  converted_range.end = gt_coords_convert_point(viewrange, node_range.end);

  return converted_range;
}
