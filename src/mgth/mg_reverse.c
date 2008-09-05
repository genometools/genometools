/*
  Copyright (c) 2007 David Schmitz-Huebsch <dschmitz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.

    modiefied for metagenomethreader;

*/

#include "mg_reverse.h"

static int mg_complement(char *reverse_char, char dna_char, Error * err)
{
  error_check(err);
  switch (dna_char)
  {
    case 'A':
      *reverse_char = 'T';
      return 0;
    case 'C':
      *reverse_char = 'G';
      return 0;
    case 'G':
      *reverse_char = 'C';
      return 0;
    case 'T':
      *reverse_char = 'A';
      return 0;
    case 'U':
      *reverse_char = 'A';
      return 0;
    case 'a':
      *reverse_char = 't';
      return 0;
    case 'c':
      *reverse_char = 'g';
      return 0;
    case 'g':
      *reverse_char = 'c';
      return 0;
    case 't':
      *reverse_char = 'a';
      return 0;
    case 'u':
      *reverse_char = 'a';
      return 0;
    case 'S':
      *reverse_char = 'S';
      return 0;
    case 's':
      *reverse_char = 's';
      return 0;
    case 'N':
      *reverse_char = 'N';
      return 0;
    case 'n':
      *reverse_char = 'n';
      return 0;
    case 'R':
      *reverse_char = 'Y';
      return 0;
    case 'r':
      *reverse_char = 'y';
      return 0;
    case 'Y':
      *reverse_char = 'R';
      return 0;
    case 'y':
      *reverse_char = 'r';
      return 0;
    case 'M':
      *reverse_char = 'K';
      return 0;
    case 'm':
      *reverse_char = 'k';
      return 0;
    case 'K':
      *reverse_char = 'M';
      return 0;
    case 'k':
      *reverse_char = 'm';
      return 0;
    case 'W':
      *reverse_char = 'W';
      return 0;
    case 'w':
      *reverse_char = 'w';
      return 0;
    case 'H':
      *reverse_char = 'D';
      return 0;
    case 'h':
      *reverse_char = 'd';
      return 0;
    case 'D':
      *reverse_char = 'H';
      return 0;
    case 'd':
      *reverse_char = 'h';
      return 0;
    case 'B':
      *reverse_char = 'V';
      return 0;
    case 'b':
      *reverse_char = 'v';
      return 0;
    case 'V':
      *reverse_char = 'B';
      return 0;
    case 'v':
      *reverse_char = 'b';
      return 0;

    default:
      error_set(err, "complement of DNA character '%c' not defined",
                dna_char);
      return -1;
  }
}

int mg_reverse_complement(char *dna_seq, unsigned long seqlen, Error * err)
{
  char *front_char,
   *back_char,
    tmp_char;
  int had_err = 0;

  error_check(err);
  assert(dna_seq);

  for (front_char = dna_seq, back_char = dna_seq + seqlen - 1;
       front_char <= back_char; front_char++, back_char--)
  {
    had_err = mg_complement(&tmp_char, *front_char, err);
    if (!had_err)
      had_err = mg_complement(front_char, *back_char, err);
    if (!had_err)
      *back_char = tmp_char;
    if (had_err)
      break;
  }
  return had_err;
}
