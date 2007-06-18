/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdbool.h>

bool islittleendian(void)
{
  int x = 1;

  if (*(char *) &x == 1)
  {
    return true;
  } else
  {
    return false;
  }
}
