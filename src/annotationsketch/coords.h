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

#ifndef COORDS_H
#define COORDS_H

#include "annotationsketch/drawing_range.h"
#include "core/range.h"

/* Returns width-independent coordinates (x,y) with 0 <= x,y <= 1
   for the given 1D coordinate. */
double         gt_coords_convert_point(GtRange viewrange, long pos);
/* Returns width-independent coordinates (x,y) with 0 <= x,y <= 1
   for the given range in the given layout. Also sets clipping info. */
GtDrawingRange gt_coords_calc_generic_range(GtRange node_range,
                                            GtRange viewrange);

#endif
