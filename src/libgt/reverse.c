/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "error.h"
#include "reverse.h"
#include "undef.h"

static char complement(char dna_char)
{
  switch (dna_char) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'a': return 't';
    case 'c': return 'g';
    case 'g': return 'c';
    case 't': return 'a';
    error("complement of DNA character '%c' not defined", dna_char);
  }
  return UNDEFCHAR; /* shut-up compiler */
}

void reverse_complement(char *dna_seq, unsigned long seqlen)
{
  char *front_char, *back_char, tmp_char = 0;
  assert(dna_seq);
  for (front_char = dna_seq, back_char = dna_seq + seqlen - 1;
      front_char <= back_char;
      front_char++, back_char--) {
    tmp_char   = complement(*front_char);
    *front_char = complement(*back_char);
    *back_char = tmp_char;
  }
}
