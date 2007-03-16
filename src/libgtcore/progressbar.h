/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

/* Start showing a progress bar on stdout if possible. */
void progressbar_start(const unsigned long *current_computation,
                       unsigned long number_of_computations);

/* Stop showing a progress bar. */
void progressbar_stop(void);

#if 0
  a typical use of the progressbar:

  i = 0;
  progressbar_start(&i, number_of_computations);
  for (; i < number_of_computations; i++) {
    /* perform the ith computation */
  }
  progressbar_stop();

#endif

#endif
