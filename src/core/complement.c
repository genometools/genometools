/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/complement.h"

int gt_complement(char *reverse_char, char dna_char, GtError *e)
{
  gt_error_check(e);
  switch (dna_char) {
    case 'A': *reverse_char = 'T'; return 0;
    case 'C': *reverse_char = 'G'; return 0;
    case 'G': *reverse_char = 'C'; return 0;
    case 'T': *reverse_char = 'A'; return 0;
    case 'N': *reverse_char = 'N'; return 0;
    case 'a': *reverse_char = 't'; return 0;
    case 'c': *reverse_char = 'g'; return 0;
    case 'g': *reverse_char = 'c'; return 0;
    case 't': *reverse_char = 'a'; return 0;
    case 'n': *reverse_char = 'n'; return 0;
    default:
      gt_error_set(e, "complement of DNA character '%c' not defined", dna_char);
      return -1;
  }
}
