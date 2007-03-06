/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

/* This implementation was heavily inspired by code from Stefan Kurtz. */

#include <assert.h>
#include <libgt/codon.h>

#define T_CODE 0
#define C_CODE 1
#define A_CODE 2
#define G_CODE 3

char codon2amino(char n0, char n1, char n2)
{
  static char aminos[] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRR"
                         "IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  unsigned int code = 0;

  switch (n0) {
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
    default: assert(0); /* XXX */
  }

  switch (n1) {
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
      code += T_CODE << 2;
      break;
    default: assert(0); /* XXX */
  }

  switch (n2) {
    case 'A':
    case 'a':
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
      code += T_CODE;
      break;
    default: assert(0); /* XXX */
  }

  return aminos[code];
}
