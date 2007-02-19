/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ARRAY2DIM_H
#define ARRAY2DIM_H

#include "xansi.h"

#define array2dim_malloc(ARRAY2DIM, ROWS, COLUMNS, TYPE)             \
        {                                                            \
          unsigned long i;                                           \
	  ARRAY2DIM = xmalloc(sizeof (TYPE*) * (ROWS));               \
	  ARRAY2DIM[0] = xmalloc(sizeof (TYPE) * (ROWS) * (COLUMNS)); \
	  for (i = 1; i < (ROWS); i++)                               \
            ARRAY2DIM[i] = ARRAY2DIM[i-1] + (COLUMNS);               \
        }

#define array2dim_calloc(ARRAY2DIM, ROWS, COLUMNS, TYPE)             \
        {                                                            \
          unsigned long i;                                           \
	  ARRAY2DIM = xmalloc(sizeof (TYPE*) * (ROWS));               \
	  ARRAY2DIM[0] = xcalloc((ROWS) * (COLUMNS), sizeof (TYPE));  \
	  for (i = 1; i < (ROWS); i++)                               \
            ARRAY2DIM[i] = ARRAY2DIM[i-1] + (COLUMNS);               \
        }

#define array2dim_free(ARRAY2DIM)                                    \
        free(ARRAY2DIM[0]);                                          \
        free(ARRAY2DIM);

#if 0
  example usage:

  double **a2dim;

  ARRAY2DIM_MALLOC(a2dim, 10, 20, double); /* create a 10 * 20 double array */
  /* ... (use array a2dim in conventional way via a2dim[row][column]) */
  ARRAY2DIM_FREE(a2dim);
#endif

#endif
