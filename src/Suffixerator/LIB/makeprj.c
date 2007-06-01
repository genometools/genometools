#include <stdbool.h>
#include <ctype.h>
#include "libgtcore/env.h"
#include "types.h"
#include "genstream.h"
#include "inputsymbol.h"
#include "chardef.h"
#include "spacedef.h"
#include "arraydef.h"
#include "fbs-def.h"
#include "dist-if.h"

#include "distcalc.pr"
#include "genericstream.pr"
#include "fbsadv.pr"

#include "readnextUchar.gen"

static Uint currentrangevalue(Uint i,Uint distvalue)
{
  if (i <= UCHAR_MAX)
  {
    return distvalue;
  }
  if (i % UCHAR_MAX == 0)
  {
    return i/UCHAR_MAX * distvalue;
  }
  return (UintConst(1) + i/UCHAR_MAX) * distvalue;
}

static int updatesumranges(void *key, void *value, void *data,
                           /*@unused@*/ Env *env)
{
  Uint keyvalue, distvalue;
  Uint *specialrangesptr = (Uint *) data;

  keyvalue = (Uint) key;
  distvalue = *((Uint *) value);
  (*specialrangesptr) += currentrangevalue(keyvalue,distvalue);
  printf("specialranges of length %lu: %lu\n",
         (Showuint) keyvalue,(Showuint) distvalue);
  return 0;
}

int scanfastasequence(
        Uint *numofsequences,
        Uint64 *totallength,
        Specialcharinfo *specialcharinfo,
        const char **filenametab,
        unsigned int numoffiles,
        PairUint **filelengthtab,
        const Uchar *symbolmap,
        Env *env)
{
  Fastabufferstate fbs;
  Uchar charcode;
  Uint64 pos;
  int retval;
  bool specialprefix = true;
  Uint lastspeciallength = 0;
  Distribution *specialrangelengths;

  *numofsequences = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;

  initfastabufferstate(&fbs,
                       filenametab,
                       numoffiles,
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
        adddistribution(specialrangelengths,lastspeciallength,env);
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
        lastspeciallength = UintConst(1);
      } else
      {
        lastspeciallength++;
      }
      if (charcode == SEPARATOR)
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
        adddistribution(specialrangelengths,lastspeciallength,env);
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
        Uint *numofsequences,
        Uint *totallength,
        PairUint **filelengthtab,
        Specialcharinfo *specialcharinfo,
        const char **filenametab,
        unsigned int numoffiles,
        const Uchar *symbolmap,
        Env *env)
{
  unsigned int filenum;
  Fgetcreturntype currentchar;
  bool indesc, firstseq = true, specialprefix = true;
  Uint linenum = UintConst(1), currentposition,
       countreadcharacters, lastspeciallength = 0;
  Genericstream inputstream;
  Uchar charcode;

  env_error_check(env);
  *numofsequences = 0;
  *totallength = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->specialranges = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  ALLOCASSIGNSPACE(*filelengthtab,NULL,PairUint,(Uint) numoffiles);
  for (filenum = 0; filenum < numoffiles; filenum++)
  {
    opengenericstream(&inputstream,filenametab[filenum]);
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
                lastspeciallength = UintConst(1);
              }
	      currentposition++;
	    }
	    indesc = true;
	  } else
	  {
	    charcode = symbolmap[(Uint) currentchar];
	    if (charcode == (Uchar) UNDEFCHAR)
	    {
              env_error_set(env,"illegal character '%c': file \"%s\", line %lu",
			    currentchar,
			    filenametab[filenum],
			    (Showuint) linenum);
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
                lastspeciallength = UintConst(1);
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
    closegenericstream(&inputstream,filenametab[filenum]);
  }
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  if (firstseq)
  {
    env_error_set(env,"no sequences in multiple fasta file(s) %s ...",
                  filenametab[0]);
    return -2;
  }
  return 0;
}
