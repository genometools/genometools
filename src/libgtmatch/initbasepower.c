#include <limits.h>
#include "spacedef.h"

unsigned int *initbasepower(unsigned int base,unsigned int len)
{
  unsigned int thepower = 1U, i, minfailure, *basepower;

  ALLOCASSIGNSPACE(basepower,NULL,unsigned int,len+1);
  minfailure = UINT_MAX/base;
  for (i=0; /* Nothing */; i++)
  {
    basepower[i] = thepower;
    if (i == len)
    {
      break;
    }
    if (thepower >= minfailure)
    {
      FREESPACE(basepower);
      fprintf(stderr,"overflow when computing %u * %u",thepower,base);
      exit(EXIT_FAILURE); /* programming error */
    }
    thepower *= base;
  }
  return basepower;
}
