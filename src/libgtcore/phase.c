/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include <libgtcore/phase.h>

Phase phase_get(char phase_char)
{
  assert(strchr(PHASECHARS, phase_char));
  switch (phase_char) {
    case '0': return PHASE_ZERO;
    case '1': return PHASE_ONE;
    case '2': return PHASE_TWO;
    default : assert(0); /*@fallthrough@*/
    case '.': return PHASE_UNDEFINED;
  }
}
