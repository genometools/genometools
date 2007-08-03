/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <ctype.h>
#include <stdbool.h>
#include <inttypes.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "fbs-def.h"
#include "stamp.h"

#include "fbsadv.pr"

#include "readnextUchar.gen"

int guessifproteinsequencestream(const StrArray *filenametab,Env *env)
{
  uint32_t countnonbases = 0,
           currentposition;
  Uchar currentchar;
  Fastabufferstate fbs;
  int retval;

  env_error_check(env);
  initformatbufferstate(&fbs,
                        filenametab,
                        NULL,
                        false,
                        NULL,
                        NULL,
                        env);
  for (currentposition = 0; currentposition < (uint32_t) 1000; 
       currentposition++)
  {
    retval = readnextUchar(&currentchar,&fbs,env);
    if (retval < 0)
    {
      return -1;
    }
    if (retval == 0)
    {
      break;
    }
    switch (currentchar)
    {
      case 'L':
      case 'I':
      case 'F':
      case 'E':
      case 'Q':
      case 'P':
      case 'X':
      case 'Z': countnonbases++;
                break;
      default:  break;
    }
  }
  if (countnonbases > 0 && countnonbases >= currentposition/10)
  {
    return 1;
  }
  return 0;
}
