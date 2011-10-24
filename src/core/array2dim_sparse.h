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

#ifndef ARRAY2DIM_SPARSE_H
#define ARRAY2DIM_SPARSE_H

#include "core/array2dim_api.h"

typedef struct {
  unsigned long offset,
                length;
} GtRowInfo;

/* Allocates a new 2-dimensional sparse array with the given number of <ROWS>
   and a total <SIZE>. Each row starts at the corresponding offset given in
   <ROWINFO> and has the corresponding length. It assigns a pointer to the newly
   allocated space to <ARRAY2DIM>. The size of each element is determined
   automatically from the type of the <ARRAY2DIM> pointer. */
#define gt_array2dim_sparse_calloc(ARRAY2DIM, ROWS, SIZE, ROWINFO)             \
        {                                                                      \
          GtRowInfo *gt_ri = (ROWINFO);                                        \
          unsigned long gt_a2d_i;                                              \
          gt_assert(gt_ri[0].offset == 0);                                     \
          ARRAY2DIM = gt_malloc(sizeof *ARRAY2DIM * (ROWS));                   \
          (ARRAY2DIM)[0] = gt_calloc((SIZE), sizeof **ARRAY2DIM);              \
          for (gt_a2d_i = 1UL; gt_a2d_i < (ROWS); gt_a2d_i++)                  \
            (ARRAY2DIM)[gt_a2d_i] = (ARRAY2DIM)[gt_a2d_i-1]                    \
                                  + gt_ri[gt_a2d_i-1].offset                   \
                                  + gt_ri[gt_a2d_i-1].length                   \
                                  - gt_ri[gt_a2d_i].offset;                    \
        }

int     gt_array2dim_sparse_example(GtError*);

#endif
