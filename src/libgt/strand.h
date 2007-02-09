/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STRAND_H
#define STRAND_H

typedef enum {
  STRAND_FORWARD, /* '+' */
  STRAND_REVERSE, /* '-' */
  STRAND_BOTH,    /* '.' */
  STRAND_UNKNOWN  /* '?' */
} Strand;

#define STRANDCHARS "+-.?"

/* returns NUM_OF_STRAND_TYPES if strand_char is not a valid one */
Strand strand_get(char strand_char);
Strand strand_join(Strand, Strand);

#endif
