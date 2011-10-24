/*
  Copyright (c) 2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/array2dim_sparse.h"
#include "core/unused_api.h"

/* allocate and fill the following sparse matrix:

              11111
    012345678901234
   +---------------+
  0|000000         |
  1| 111111        |
  2|  222222       |
  3|   333333      |
  4|    444444     |
  5|     555555    |
  6|      666666   |
  7|       777777  |
  8|        888888 |
  9|         999999|
   +---------------+
*/
int gt_array2dim_sparse_example(GT_UNUSED GtError *err)
{
  int **a2dim, i, j;
  GtRowInfo ri[10];
  gt_error_check(err);

  /* initialize row info */
  for (i = 0; i < 10; i++) {
    ri[i].offset = i;
    ri[i].length = 6;
  }

  /* create sparse matrix */
  gt_array2dim_sparse_calloc(a2dim, 10, 60, ri);

  /* fill matrix */
  for (i = 0; i < 10; i++) {
    for (j = 0; j < 6; j++)
      a2dim[i][i+j] = i;
  }

  gt_assert(a2dim[0][0] == 0);
  gt_assert(a2dim[5][5] == 5);
  gt_assert(a2dim[9][14] == 9);

  /* free */
  gt_array2dim_delete(a2dim);

  return 0;
}
