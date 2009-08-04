/*
  Copyright (c) 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/array3dim.h"
#include "core/unused_api.h"

/* example usage of the array3dim macros */
int gt_array3dim_example(GT_UNUSED GtError *err)
{
  double ***a3dim;
  int x, y, z;
  gt_error_check(err);

  /* create a 10 x 20 x 30 double array */
  gt_array3dim_malloc(a3dim, 10, 20, 30);

  /* ... (use array a3dim in conventional way via a3dim[x][y][z]) */
  for (x = 1; x < 10; x++) {
    for (y = 1; y < 20; y++) {
      for (z = 1; z < 30; z++)
        a3dim[x][y][z] = x + y + z;
    }
  }

  /* free */
  gt_array3dim_delete(a3dim);

  return 0;
}
