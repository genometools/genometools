/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "intcode-def.h"

void kmercode2string(char *buffer,
                     Codetype code,
                     uint32_t numofchars,
                     uint32_t kmersize,
                     const char *characters)
{
  int i;
  uint32_t cc; 
  Codetype tmpcode = code;

  buffer[kmersize] = '\0';
  for (i=(int) (kmersize-1); i>=0; i--)
  {
    cc = tmpcode % numofchars;
    buffer[i] = (char) characters[cc];
    tmpcode = (tmpcode - cc) / numofchars;
  }
}
