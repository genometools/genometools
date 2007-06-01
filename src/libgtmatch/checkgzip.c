/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <string.h>

#define GZIPSUFFIX        ".gz"
#define GZIPSUFFIXLENGTH  (sizeof (GZIPSUFFIX)-1)

unsigned char checkgzipsuffix(const char *filename)
{
  size_t filenamelength = strlen (filename);

  if (filenamelength < GZIPSUFFIXLENGTH ||
      strcmp (filename + filenamelength - GZIPSUFFIXLENGTH, GZIPSUFFIX) != 0)
  {
    return 0;
  }
  return (unsigned char) 1;
}
