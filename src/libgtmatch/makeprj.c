/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdbool.h>
#include <ctype.h>
#include "libgtcore/env.h"
#include "libgtcore/strarray.h"
#include "types.h"
#include "genstream.h"
#include "inputsymbol.h"
#include "chardef.h"
#include "spacedef.h"
#include "arraydef.h"
#include "fbs-def.h"
#include "dist-if.h"
#include "safecast-gen.h"
#include "stamp.h"

#include "distcalc.pr"
#include "genericstream.pr"
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
  printf("specialranges of length %lu: %lu\n",keyvalue,distvalue);
  return 0;
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

int scanfastasequence(
        unsigned long *numofsequences,
        Seqpos *totallength,
        Specialcharinfo *specialcharinfo,
        const StrArray *filenametab,
        PairSeqpos **filelengthtab,
        const Uchar *symbolmap,
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

  *numofsequences = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;

  initfastabufferstate(&fbs,
                       filenametab,
                       symbolmap,
                       filelengthtab,
                       env);
  specialrangelengths = initdistribution(env);
  for (pos = 0; /* Nothing */; pos++)
  {
    retval = readnextUchar(&charcode,&fbs,env);
    if (retval < 0)
    {
      return -1;
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
  specialcharinfo->specialranges = 0;
  (void) foreachdistributionvalue(specialrangelengths,updatesumranges,
                                  &specialcharinfo->specialranges,env);
  freedistribution(&specialrangelengths,env);
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  (*numofsequences)++;
  *totallength = pos;
  return 0;
}

/* the following function is obsolete */

int scanfastasequence2(
        unsigned long *numofsequences,
        Seqpos *totallength,
        PairSeqpos **filelengthtab,
        Specialcharinfo *specialcharinfo,
        const StrArray *filenametab,
        const Uchar *symbolmap,
        Env *env)
{
  unsigned long filenum;
  uint32_t linenum = (uint32_t) 1;
  Fgetcreturntype currentchar;
  bool indesc, firstseq = true, specialprefix = true;
  Seqpos currentposition,
       countreadcharacters,
       lastspeciallength = 0;
  Genericstream inputstream;
  Uchar charcode;

  env_error_check(env);
  *numofsequences = 0;
  *totallength = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->specialranges = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  ALLOCASSIGNSPACE(*filelengthtab,NULL,PairSeqpos,
                   strarray_size(filenametab));
  for (filenum = 0; filenum < strarray_size(filenametab); filenum++)
  {
    opengenericstream(&inputstream,strarray_get(filenametab,filenum));
    indesc = false;
    currentposition = 0;
    countreadcharacters = 0;
    for (countreadcharacters = 0; /* Nothing */; countreadcharacters++)
    {
      if (inputstream.isgzippedstream)
      {
	currentchar = gzgetc(inputstream.stream.gzippedstream);
      } else
      {
	currentchar = fgetc(inputstream.stream.fopenstream);
      }
      if (currentchar == EOF)
      {
	break;
      }
      if (indesc)
      {
	if (currentchar == NEWLINESYMBOL)
	{
	  linenum++;
	  indesc = false;
	}
      } else
      {
	if (!isspace((Ctypeargumenttype) currentchar))
	{
	  if (currentchar == FASTASEPARATOR)
	  {
            (*numofsequences)++;
	    if (firstseq)
	    {
	      firstseq = false;
	    } else
	    {
              specialcharinfo->specialcharacters++;
              if (lastspeciallength == 0)
              {
                specialcharinfo->specialranges++;
                lastspeciallength = (Seqpos) 1;
              }
	      currentposition++;
	    }
	    indesc = true;
	  } else
	  {
	    charcode = symbolmap[(uint32_t) currentchar];
	    if (charcode == (Uchar) UNDEFCHAR)
	    {
              env_error_set(env,"illegal character '%c': file \"%s\", line %u",
			    currentchar,
			    strarray_get(filenametab,filenum),
			    (unsigned int) linenum);
	      return -1;
	    }
	    currentposition++;
	    if (ISSPECIAL(charcode))
	    {
              if (specialprefix)
              {
                specialcharinfo->lengthofspecialprefix++;
              }
              specialcharinfo->specialcharacters++;
              if (lastspeciallength == 0)
              {
                specialcharinfo->specialranges++;
                lastspeciallength = (Seqpos) 1;
              } else
              {
                lastspeciallength++;
              }
	    } else
            {
              specialprefix = false;
              if (lastspeciallength > 0)
              {
                lastspeciallength = 0;
              }
            }
	  }
	}
      }
    }
    (*filelengthtab)[filenum].uint0 = countreadcharacters;
    (*filelengthtab)[filenum].uint1 = currentposition;
    *totallength += currentposition;
    closegenericstream(&inputstream,strarray_get(filenametab,filenum));
  }
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  if (firstseq)
  {
    env_error_set(env,"no sequences in multiple fasta file(s) %s ...",
                  strarray_get(filenametab,0));
    return -2;
  }
  return 0;
}
