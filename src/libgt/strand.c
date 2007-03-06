/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include <libgt/strand.h>

Strand strand_get(char strand_char)
{
  assert(strchr(STRANDCHARS, strand_char));
  switch (strand_char) {
    case '+': return STRAND_FORWARD;
    case '-': return STRAND_REVERSE;
    case '.': return STRAND_BOTH;
    default : assert(0); /*@fallthrough@*/
    case '?': return STRAND_UNKNOWN;
  }
}

Strand strand_join(Strand strand_a, Strand strand_b)
{
  switch (strand_b) {
    case STRAND_FORWARD:
      assert(strand_a != STRAND_REVERSE);
      return STRAND_FORWARD;
    case STRAND_REVERSE:
      assert(strand_a != STRAND_FORWARD);
      return STRAND_REVERSE;
    case STRAND_BOTH:
    case STRAND_UNKNOWN:
      /* strand_a == STRAND_FORWARD -> stays the same */
      /* strand_a == STRAND_REVERSE -> stays the same */
      /* strand_a == STRAND_UNKNOWN -> stays the same */
      if (strand_a == STRAND_BOTH)
        return STRAND_UNKNOWN;
  }
  return strand_a;
}
