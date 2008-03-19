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
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

/* This implementation was heavily inspired by code from Stefan Kurtz. */

/* erweitert fuer MG um den Fall unbekannter AS X und den Faellen
   der Degeneration an der dritten Stelle (R/Y/N) */

/*     modified for metagenomethreader; */

#include "mg_codon.h"

#define T_CODE 0
#define C_CODE 1
#define A_CODE 2
#define G_CODE 3
#define X_CODE 4

char mg_codon2amino(char n0, char n1, char n2)
{
  static char aminos[] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGGX";
  unsigned int code = 0;

  switch (n0)
  {
    case 'A':
    case 'a':
      code = A_CODE << 4;
      break;
    case 'C':
    case 'c':
      code = C_CODE << 4;
      break;
    case 'G':
    case 'g':
      code = G_CODE << 4;
      break;
    case 'T':
    case 't':
    case 'U':
    case 'u':
      code = T_CODE << 4;
      break;
    default:
      code = G_CODE << 4;
  }

  switch (n1)
  {
    case 'A':
    case 'a':
      code += A_CODE << 2;
      break;
    case 'C':
    case 'c':
      code += C_CODE << 2;
      break;
    case 'G':
    case 'g':
      code += G_CODE << 2;
      break;
    case 'T':
    case 't':
    case 'U':
    case 'u':
      /* code += T_CODE << 2; */
      break;
    default:
      code = G_CODE << 2;
  }

  switch (n2)
  {
    case 'A':
    case 'a':
    case 'R':
    case 'r':
      code += A_CODE;
      break;
    case 'C':
    case 'c':
      code += C_CODE;
      break;
    case 'G':
    case 'g':
      code += G_CODE;
      break;
    case 'T':
    case 't':
    case 'U':
    case 'u':
    case 'N':
    case 'n':
    case 'Y':
    case 'y':
      code += T_CODE;
      break;
    default:
      code += X_CODE;
  }

  return aminos[code];
}
