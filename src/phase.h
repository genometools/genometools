/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef PHASE_H
#define PHASE_H

typedef enum {
  PHASE_ZERO,     /* '0' */
  PHASE_ONE,      /* '1' */
  PHASE_TWO,      /* '2' */
  PHASE_UNDEFINED /* '.' */
} Phase;

#define PHASECHARS "012."

/* returns NUM_OF_PHASES if phase_char is not a valid one */
Phase phase_get(char phase_char);

#endif
