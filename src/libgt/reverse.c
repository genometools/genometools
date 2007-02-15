/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "error.h"
#include "reverse.h"
#include "undef.h"

static int complement(char *reverse_char, char dna_char, Error *err)
{
  error_check(err);
  switch (dna_char) {
    case 'A': *reverse_char = 'T'; return 0;
    case 'C': *reverse_char = 'G'; return 0;
    case 'G': *reverse_char = 'C'; return 0;
    case 'T': *reverse_char = 'A'; return 0;
    case 'a': *reverse_char = 't'; return 0;
    case 'c': *reverse_char = 'g'; return 0;
    case 'g': *reverse_char = 'c'; return 0;
    case 't': *reverse_char = 'a'; return 0;
    default:
      error_set(err, "complement of DNA character '%c' not defined", dna_char);
      return -1;
  }
}

int reverse_complement(char *dna_seq, unsigned long seqlen, Error *err)
{
  char *front_char, *back_char, tmp_char;
  int has_err = 0;
  error_check(err);
  assert(dna_seq);
  for (front_char = dna_seq, back_char = dna_seq + seqlen - 1;
       front_char <= back_char;
       front_char++, back_char--) {
    has_err = complement(&tmp_char, *front_char, err);
    if (!has_err)
      has_err = complement(front_char, *back_char, err);
    if (!has_err)
      *back_char = tmp_char;
    if (has_err)
      break;
  }
  return has_err;
}
