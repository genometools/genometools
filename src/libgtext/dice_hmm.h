/*
  Copyright (c) 2005-2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef DICE_HMM

#include <libgtext/hmm.h>

typedef enum {
  DICE_FAIR,
  DICE_LOADED,
  DICE_NUM_OF_STATES
} Dice_states;

typedef enum {
  ONE,
  TWO,
  THREE,
  FOUR,
  FIVE,
  SIX,
  DICE_NUM_OF_SYMBOLS
} Dice_emissions;

HMM*   dice_hmm_loaded(Env*);
HMM*   dice_hmm_fair(Env*);
Alpha* dice_hmm_alpha(Env*);

#endif
