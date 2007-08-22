/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdbool.h>
#include <ctype.h>
#include "libgtcore/env.h"
#include "libgtcore/strarray.h"
#include "seqpos-def.h"
#include "chardef.h"
#include "spacedef.h"
#include "arraydef.h"
#include "fbs-def.h"
#include "dist-if.h"
#include "safecast-gen.h"
#include "stamp.h"

#include "distcalc.pr"
#include "fbsadv.pr"

#include "readnextUchar.gen"

static unsigned long currentrangevalue(unsigned long i,unsigned long distvalue)
{
  if (i <= UCHAR_MAX)
  {
    return distvalue;
  }
  if (i % UCHAR_MAX == 0)
  {
    return i/UCHAR_MAX * distvalue;
  }
  return (((unsigned long) 1) + i/UCHAR_MAX) * distvalue;
}

static int updatesumranges(void *key, void *value, void *data,
                           /*@unused@*/ Env *env)
{
  unsigned long keyvalue,
                distvalue,
                *specialrangesptr = (unsigned long *) data;

  keyvalue = (unsigned long) key;
  distvalue = *((unsigned long *) value);
  (*specialrangesptr) += currentrangevalue(keyvalue,distvalue);
  printf("# specialranges of length %lu: %lu\n",keyvalue,distvalue);
  return 0;
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

int scanfastasequence(
        unsigned long *numofsequences,
        Seqpos *totallength,
        Specialcharinfo *specialcharinfo,
        const StrArray *filenametab,
        Filelengthvalues **filelengthtab,
        const Uchar *symbolmap,
        bool plainformat,
        Env *env)
{
  Fastabufferstate fbs;
  Uchar charcode;
  Seqpos pos;
  int retval;
  bool specialprefix = true;
  Seqpos lastspeciallength = 0;
  Distribution *specialrangelengths;
  unsigned long idx;
  bool haserr = false;

  env_error_check(env);
  *numofsequences = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;

  initformatbufferstate(&fbs,
                        filenametab,
                        symbolmap,
                        plainformat,
                        filelengthtab,
                        NULL,
                        env);
  specialrangelengths = initdistribution(env);
  for (pos = 0; /* Nothing */; pos++)
  {
    retval = readnextUchar(&charcode,&fbs,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      if (lastspeciallength > 0)
      {
        idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
        adddistribution(specialrangelengths,idx,env);
      }
      break;
    }
    if (ISSPECIAL(charcode))
    {
      if (specialprefix)
      {
        specialcharinfo->lengthofspecialprefix++;
      }
      specialcharinfo->specialcharacters++;
      if (lastspeciallength == 0)
      {
        lastspeciallength = (Seqpos) 1;
      } else
      {
        lastspeciallength++;
      }
      if (charcode == (Uchar) SEPARATOR)
      {
        (*numofsequences)++;
      }
    } else
    {
      if (specialprefix)
      {
        specialprefix = false;
      }
      if (lastspeciallength > 0)
      {
        idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
        adddistribution(specialrangelengths,idx,env);
        lastspeciallength = 0;
      }
    }
  }
  if (!haserr)
  {
    specialcharinfo->specialranges = 0;
    (void) foreachdistributionvalue(specialrangelengths,updatesumranges,
                                    &specialcharinfo->specialranges,env);
    specialcharinfo->lengthofspecialsuffix = lastspeciallength;
    (*numofsequences)++;
    *totallength = pos;
  }
  freedistribution(&specialrangelengths,env);
  return haserr ? -1 : 0;
}
