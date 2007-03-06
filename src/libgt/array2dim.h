/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ARRAY2DIM_H
#define ARRAY2DIM_H

#include <libgt/xansi.h>

#define array2dim_malloc(ARRAY2DIM, ROWS, COLUMNS, TYPE, ENV)               \
        {                                                                   \
          unsigned long i;                                                  \
	  ARRAY2DIM = env_ma_malloc(ENV, sizeof (TYPE*) * (ROWS));          \
	  (ARRAY2DIM)[0] = env_ma_malloc(ENV,                               \
                                       sizeof (TYPE) * (ROWS) * (COLUMNS)); \
	  for (i = 1; i < (ROWS); i++)                                      \
            (ARRAY2DIM)[i] = (ARRAY2DIM)[i-1] + (COLUMNS);                  \
        }

#define array2dim_calloc(ARRAY2DIM, ROWS, COLUMNS, TYPE, ENV)               \
        {                                                                   \
          unsigned long i;                                                  \
	  ARRAY2DIM = env_ma_malloc(ENV, sizeof (TYPE*) * (ROWS));          \
	  (ARRAY2DIM)[0] = env_ma_calloc(ENV, (ROWS) * (COLUMNS),           \
                                       sizeof (TYPE));                      \
	  for (i = 1; i < (ROWS); i++)                                      \
            (ARRAY2DIM)[i] = (ARRAY2DIM)[i-1] + (COLUMNS);                  \
        }

#define array2dim_delete(ARRAY2DIM, ENV)                                    \
        env_ma_free((ARRAY2DIM)[0], ENV);                                   \
        env_ma_free(ARRAY2DIM, ENV);

#if 0
  example usage:

  double **a2dim;

  /* create a 10 * 20 double array */
  array2dim_malloc(a2dim, 10, 20, double, env);
  /* ... (use array a2dim in conventional way via a2dim[row][column]) */
  array2dim_delete(a2dim, env);
#endif

#endif
