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

#ifndef ARRAY3DIM_H
#define ARRAY3DIM_H

#include "core/error_api.h"
#include "core/ma_api.h"

/* Array3dim module */

/* Allocates a new 3-dimensional array with dimensions <X_SIZE> x <Y_SIZE> x
   <Z_SIZE> and assigns a pointer to the newly allocated space to <ARRAY3DIM>.
   The size of each element is determined automatically from the type of the
   <ARRAY3DIM> pointer. */
#define gt_array3dim_malloc(ARRAY3DIM, X_SIZE, Y_SIZE, Z_SIZE)                 \
        {                                                                      \
          unsigned long gt_a3d_x, gt_a3d_y;                                    \
          ARRAY3DIM = gt_malloc(sizeof *ARRAY3DIM * (X_SIZE));                 \
          (ARRAY3DIM)[0] = gt_malloc(sizeof **ARRAY3DIM * (X_SIZE) * (Y_SIZE));\
          for (gt_a3d_x = 1UL; gt_a3d_x < (X_SIZE); gt_a3d_x++)                \
            (ARRAY3DIM)[gt_a3d_x] = (ARRAY3DIM)[gt_a3d_x-1] + (Y_SIZE);        \
          (ARRAY3DIM)[0][0] = gt_malloc(sizeof ***ARRAY3DIM *                  \
                                        (X_SIZE) * (Y_SIZE) * (Z_SIZE));       \
          for (gt_a3d_y = 1UL; gt_a3d_y < (Y_SIZE); gt_a3d_y++) {              \
            (ARRAY3DIM)[0][gt_a3d_y] = ARRAY3DIM[0][gt_a3d_y-1] + (Z_SIZE);    \
          }                                                                    \
          for (gt_a3d_x = 1UL; gt_a3d_x < (X_SIZE); gt_a3d_x++) {              \
            (ARRAY3DIM)[gt_a3d_x][0] =                                         \
              (ARRAY3DIM)[gt_a3d_x-1][(Y_SIZE)-1] + (Z_SIZE);                  \
            for (gt_a3d_y = 1UL; gt_a3d_y < (Y_SIZE); gt_a3d_y++) {            \
              (ARRAY3DIM)[gt_a3d_x][gt_a3d_y] =                                \
                (ARRAY3DIM)[gt_a3d_x][gt_a3d_y-1] + (Z_SIZE);                  \
            }                                                                  \
          }                                                                    \
        }

/* Allocates a new 3-dimensional array with dimensions <X_SIZE> x <Y_SIZE> and
   assigns a pointer to the newly allocated space to <ARRAY3DIM>.
   The allocated space is initialized to be filled with zeroes.
   The size of each element is determined automatically from the type of the
   <ARRAY3DIM> pointer. */
#define gt_array3dim_calloc(ARRAY3DIM, X_SIZE, Y_SIZE, Z_SIZE)                 \
        {                                                                      \
          unsigned long gt_a3d_x, gt_a3d_y;                                    \
          ARRAY3DIM = gt_malloc(sizeof *ARRAY3DIM * (X_SIZE));                 \
          (ARRAY3DIM)[0] = gt_malloc(sizeof **ARRAY3DIM * (X_SIZE) * (Y_SIZE));\
          for (gt_a3d_x = 1UL; gt_a3d_x < (X_SIZE); gt_a3d_x++)                \
            (ARRAY3DIM)[gt_a3d_x] = (ARRAY3DIM)[gt_a3d_x-1] + (Y_SIZE);        \
          (ARRAY3DIM)[0][0] = gt_calloc(sizeof ***ARRAY3DIM *                  \
                                        (X_SIZE) * (Y_SIZE) * (Z_SIZE));       \
          for (gt_a3d_y = 1UL; gt_a3d_y < (Y_SIZE); gt_a3d_y++) {              \
            (ARRAY3DIM)[0][gt_a3d_y] = ARRAY3DIM[0][gt_a3d_y-1] + (Z_SIZE);    \
          }                                                                    \
          for (gt_a3d_x = 1UL; gt_a3d_x < (X_SIZE); gt_a3d_x++) {              \
            (ARRAY3DIM)[gt_a3d_x][0] =                                         \
              (ARRAY3DIM)[gt_a3d_x-1][(Y_SIZE)-1] + (Z_SIZE);                  \
            for (gt_a3d_y = 1UL; gt_a3d_y < (Y_SIZE); gt_a3d_y++) {            \
              (ARRAY3DIM)[gt_a3d_x][gt_a3d_y] =                                \
                (ARRAY3DIM)[gt_a3d_x][gt_a3d_y-1] + (Z_SIZE);                  \
            }                                                                  \
          }                                                                    \
        }

/* An example for usage of the <Array3dim> module. */
int     gt_array3dim_example(GtError*);

/* Frees the space allocated for the 3-dimensional array pointed to by
   <ARRAY3DIM>. */
#define gt_array3dim_delete(ARRAY3DIM) \
        gt_free((ARRAY3DIM)[0][0]);    \
        gt_free((ARRAY3DIM)[0]);       \
        gt_free(ARRAY3DIM);

#endif
