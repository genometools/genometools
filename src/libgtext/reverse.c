/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <assert.h>
#include "libgtext/reverse.h"

static int complement(char *reverse_char, char dna_char, Error *e)
{
  error_check(e);
  switch (dna_char) {
    case 'A': *reverse_char = 'T'; return 0;
    case 'C': *reverse_char = 'G'; return 0;
    case 'G': *reverse_char = 'C'; return 0;
    case 'T': *reverse_char = 'A'; return 0;
    case 'a': *reverse_char = 't'; return 0;
    case 'c': *reverse_char = 'g'; return 0;
    case 'g': *reverse_char = 'c'; return 0;
    case 't': *reverse_char = 'a'; return 0;
    case 'n': *reverse_char = 'n'; return 0;
    default:
      error_set(e, "complement of DNA character '%c' not defined", dna_char);
      return -1;
  }
}

int reverse_complement(char *dna_seq, unsigned long seqlen, Error *e)
{
  char *front_char, *back_char, tmp_char;
  int had_err = 0;
  error_check(e);
  assert(dna_seq);
  for (front_char = dna_seq, back_char = dna_seq + seqlen - 1;
       front_char <= back_char;
       front_char++, back_char--) {
    had_err = complement(&tmp_char, *front_char, e);
    if (!had_err)
      had_err = complement(front_char, *back_char, e);
    if (!had_err)
      *back_char = tmp_char;
    if (had_err)
      break;
  }
  return had_err;
}
