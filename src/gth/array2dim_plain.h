/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef ARRAY2DIM_PLAIN_H
#define ARRAY2DIM_PLAIN_H

/* Special case of the gt_array2dim_malloc() macro taken from core/arrayd2dim.h
   which uses a plain malloc(3) instead of gt_malloc() and can therefore fail
   (i.e., return NULL). */

#define gth_array2dim_plain_malloc(ARRAY2DIM, ROWS, COLUMNS)                   \
        {                                                                      \
          unsigned long gth_a2d_i;                                             \
          ARRAY2DIM = malloc(sizeof *ARRAY2DIM * (ROWS));                      \
          if (ARRAY2DIM) {                                                     \
            (ARRAY2DIM)[0] = malloc(sizeof **ARRAY2DIM * (ROWS) * (COLUMNS));  \
            if ((ARRAY2DIM)[0]) {                                              \
              for (gth_a2d_i = 1; gth_a2d_i < (ROWS); gth_a2d_i++)             \
                (ARRAY2DIM)[gth_a2d_i] = (ARRAY2DIM)[gth_a2d_i-1] + (COLUMNS); \
            }                                                                  \
            else {                                                             \
              free(ARRAY2DIM);                                                 \
              ARRAY2DIM = NULL;                                                \
            }                                                                  \
          }                                                                    \
        }

/* Special case of the gt_array2dim_calloc() macro taken from core/arrayd2dim.h
   which uses a plain calloc(3) instead of gt_calloc() and can therefore fail
   (i.e., return NULL). */

#define gth_array2dim_plain_calloc(ARRAY2DIM, ROWS, COLUMNS)                   \
        {                                                                      \
          unsigned long gth_a2d_i;                                             \
          ARRAY2DIM = calloc((ROWS), sizeof *ARRAY2DIM);                       \
          if (ARRAY2DIM) {                                                     \
            (ARRAY2DIM)[0] = calloc((ROWS) * (COLUMNS), sizeof **ARRAY2DIM);   \
            if ((ARRAY2DIM)[0]) {                                              \
              for (gth_a2d_i = 1; gth_a2d_i < (ROWS); gth_a2d_i++)             \
                (ARRAY2DIM)[gth_a2d_i] = (ARRAY2DIM)[gth_a2d_i-1] + (COLUMNS); \
            }                                                                  \
            else {                                                             \
              free(ARRAY2DIM);                                                 \
              ARRAY2DIM = NULL;                                                \
            }                                                                  \
          }                                                                    \
        }

#define gth_array2dim_plain_delete(ARRAY2DIM) \
        free((ARRAY2DIM)[0]);                 \
        free(ARRAY2DIM);

#endif
