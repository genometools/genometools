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

#include <ctype.h>
#include "core/complement.h"

int gt_complement(char *reverse_char, char dna_char, GtError *err)
{
  gt_error_check(err);
  switch (dna_char) {
    case 'A': *reverse_char = 'T'; return 0; /* adenine */
    case 'T': *reverse_char = 'A'; return 0; /* thymine */
    case 'U': *reverse_char = 'A'; return 0; /* uracil */
    case 'G': *reverse_char = 'C'; return 0; /* guanine */
    case 'C': *reverse_char = 'G'; return 0; /* cytosine */
    case 'Y': *reverse_char = 'R'; return 0; /* pyrimidine */
    case 'R': *reverse_char = 'Y'; return 0; /* purine */
    case 'S': *reverse_char = 'S'; return 0; /* strong (3 hbonds) */
    case 'W': *reverse_char = 'W'; return 0; /* weak (2 hbonds)  */
    case 'K': *reverse_char = 'M'; return 0; /* keto */
    case 'M': *reverse_char = 'K'; return 0; /* amino */
    case 'B': *reverse_char = 'V'; return 0; /* not A */
    case 'D': *reverse_char = 'H'; return 0; /* not C */
    case 'H': *reverse_char = 'D'; return 0; /* not G */
    case 'V': *reverse_char = 'B'; return 0; /* not T/U */
    case 'N': *reverse_char = 'N'; return 0; /* any */
    case 'a': *reverse_char = 't'; return 0; /* adenine */
    case 't': *reverse_char = 'a'; return 0; /* thymine */
    case 'u': *reverse_char = 'a'; return 0; /* uracil */
    case 'g': *reverse_char = 'c'; return 0; /* guanine */
    case 'c': *reverse_char = 'g'; return 0; /* cytosine */
    case 'y': *reverse_char = 'r'; return 0; /* pyrimidine */
    case 'r': *reverse_char = 'y'; return 0; /* purine */
    case 's': *reverse_char = 's'; return 0; /* strong (3 hbonds) */
    case 'w': *reverse_char = 'w'; return 0; /* weak (2 hbonds)  */
    case 'k': *reverse_char = 'm'; return 0; /* keto */
    case 'm': *reverse_char = 'k'; return 0; /* amino */
    case 'b': *reverse_char = 'v'; return 0; /* not A */
    case 'd': *reverse_char = 'h'; return 0; /* not C */
    case 'h': *reverse_char = 'd'; return 0; /* not G */
    case 'v': *reverse_char = 'b'; return 0; /* not T/U */
    case 'n': *reverse_char = 'n'; return 0; /* any */
    default:
      if (isspace(dna_char)) {
        gt_error_set(err, "complement of whitespace character not defined");
      } else {
        gt_error_set(err, "complement of DNA character '%c' not defined",
                     dna_char);
      }
      return -1;
  }
}
