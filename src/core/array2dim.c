/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/array2dim_api.h"
#include "core/unused_api.h"

/* example usage of the array2dim macros */
int gt_array2dim_example(GT_UNUSED GtError *err)
{
  double **a2dim;
  int i, j;
  gt_error_check(err);

  /* create a 10 x 20 double array */
  gt_array2dim_malloc(a2dim, 10, 20);

  /* ... (use array a2dim in conventional way via a2dim[row][column]) */
  for (i = 1; i < 10; i++) {
    for (j = 1; j < 20; j++)
      a2dim[i][j] = i + j;
  }

  /* free */
  gt_array2dim_delete(a2dim);

  return 0;
}
