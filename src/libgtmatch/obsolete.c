/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

/* This file contains obsolete code */

/* from mappedstr.c

int getfastastreamkmers2(
        const StrArray *filenametab,
        void(*processkmercode)(void *,Codetype,Seqpos,
                               const Firstspecialpos *,Env *),
        void *processkmercodeinfo,
        uint32_t numofchars,
        uint32_t kmersize,
        const Uchar *symbolmap,
        Env *env)
{
  unsigned long filenum;
  Fgetcreturntype currentchar;
  bool indesc, firstseq = true;
  uint32_t overshoot;
  uint32_t linenum = (uint32_t) 1;
  Seqpos currentposition = 0;
  Streamstate spwp;
  Genericstream inputstream;
  Uchar charcode;

  env_error_check(env);
  initmultimappower(&spwp.multimappower,numofchars,kmersize,env);
  spwp.lengthwithoutspecial = 0;
  spwp.codewithoutspecial = 0;
  spwp.kmersize = kmersize;
  spwp.numofchars = numofchars;
  spwp.windowwidth = 0;
  spwp.firstindex = 0;
  ALLOCASSIGNSPACE(spwp.cyclicwindow,NULL,Uchar,kmersize);
  specialemptyqueue(&spwp.spos,kmersize,env);
  filllargestchartable(&spwp.filltable,numofchars,kmersize,env);
  for (filenum = 0; filenum < strarray_size(filenametab); filenum++)
  {
    opengenericstream(&inputstream,strarray_get(filenametab,filenum));
    indesc = false;
    for (;;)
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
            if (firstseq)
            {
              firstseq = false;
            } else
            {
              shiftrightwithchar(processkmercode,processkmercodeinfo,
                                 &spwp,currentposition,(Uchar) SEPARATOR,env);
              currentposition++;
            }
            indesc = true;
          } else
          {
            charcode = symbolmap[(uint32_t) currentchar];
            if (charcode == (Uchar) UNDEFCHAR)
            {
              env_error_set(env,
                            "illegal character '%c': file \"%s\", line %u",
                            currentchar,
                            strarray_get(filenametab,filenum),
                            (unsigned int) linenum);
              return -1;
            }
            shiftrightwithchar(processkmercode,processkmercodeinfo,
                               &spwp,currentposition,charcode,env);
            currentposition++;
          }
        }
      }
    }
    closegenericstream(&inputstream,strarray_get(filenametab,filenum));
  }
  for (overshoot=0; overshoot<kmersize; overshoot++)
  {
    shiftrightwithchar(processkmercode,processkmercodeinfo,&spwp,
                       currentposition + overshoot,(Uchar) WILDCARD,env);
  }
  FREESPACE(spwp.cyclicwindow);
  FREESPACE(spwp.filltable);
  ARRAY2DIMFREE(spwp.multimappower);
  specialwrapqueue(&spwp.spos,env);
  if (firstseq)
  {
    env_error_set(env,"no sequences in multiple fasta file(s) %s ...",
                  strarray_get(filenametab,0));
    return -2;
  }
  return 0;
}
*/

/* from completespecials

typedef struct
{
  uint32_t prefixlength;
  Seqpos totallength,
         countspecialmaxprefixlen0;
} CountCompletespecials;

static int lengthofspecialranges(void *info,const PairSeqpos *pair,Env *env)
{
  uint32_t len = (uint32_t) (pair->uint1 - pair->uint0);
  CountCompletespecials *csp = (CountCompletespecials *) info;

  if (pair->uint0 == 0)
  {
    if (pair->uint1 == csp->totallength)
    {
      csp->countspecialmaxprefixlen0 += (len+1);
    } else
    {
      csp->countspecialmaxprefixlen0 += len;
    }
  } else
  {
    if (pair->uint1 == csp->totallength)
    {
      csp->countspecialmaxprefixlen0 += len;
    } else
    {
      if (len >= (uint32_t) 2)
      {
        csp->countspecialmaxprefixlen0 += (len - 1);
      }
    }
  }
  return 0;
}

Seqpos determinefullspecials(const Encodedsequence *encseq,
                             Seqpos totallength,
                             uint32_t prefixlength,
                             Env *env)
{
  CountCompletespecials csp;

  csp.countspecialmaxprefixlen0 = 0;
  csp.prefixlength = prefixlength;
  csp.totallength = totallength;
  (void) overallspecialranges(encseq,lengthofspecialranges,&csp,env);
  return csp.countspecialmaxprefixlen0;
}
*/

/* from compfilename.c

char *composefilenamegeneric(const char *file,
                                           int linenum,
                                           const char *filename,
                                           char sep,
                                           const char *suffix,
                                           Env *env)
{
  size_t i, lenfilename, lensuffix, totalsize;
  char *dest;

  assert(filename != NULL);
  lenfilename = strlen(filename);
  assert(suffix != NULL);
  lensuffix = strlen(suffix);
  totalsize = lenfilename+lensuffix+1+1;
  ALLOCASSIGNSPACEGENERIC(file,linenum,dest,NULL,char,totalsize);
  assert(dest != NULL);
  for (i=0; i<lenfilename; i++)
  {
    dest[i] = filename[i];
  }
  dest[lenfilename] = sep;
  for (i=0; i<lensuffix; i++)
  {
    dest[lenfilename+1+i] = suffix[i];
  }
  dest[lenfilename+lensuffix+1] = '\0';
  return dest;
}

char *composefilename(const char *file,
                                    int linenum,
                                    const char *filename,
                                    const char *suffix,
                                    Env *env)
{
  return composefilenamegeneric(file,
                                linenum,
                                filename,
                                '.',
                                suffix,
                                env);
}
*/

/* from dstrdup.c

char *dynamicstrdup(const char *file,int linenum,
                                  const char *source,Env *env)
{
  size_t sourcelength;
  char *dest;

  assert(source != NULL);
  sourcelength = strlen(source);
  ALLOCASSIGNSPACEGENERIC(file,linenum,dest,NULL,char,sourcelength+1);
  strcpy(dest,source);
  return dest;
}
*/

/* from makeprj.c

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
    for (countreadcharacters = 0; ; countreadcharacters++)
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
*/

/* from guessprot

bool guessifproteinsequence(const Uchar *input,Seqpos inputlen)
{
  Uchar current;
  const Uchar *inputptr;
  Seqpos countnonbases = 0,
         countcharacters = 0,
         readnumoffirstcharacters,
         leastnumofnonbases;
  bool indesc = false;

  if (inputlen < (Seqpos) 1000)
  {
    readnumoffirstcharacters = inputlen;
  } else
  {
    readnumoffirstcharacters = (Seqpos) 1000;
  }
  leastnumofnonbases = readnumoffirstcharacters/10;
  for (inputptr = input; countnonbases < leastnumofnonbases &&
                        countcharacters < readnumoffirstcharacters &&
                        inputptr < input + inputlen; inputptr++)
  {
    current = *inputptr;
    if (indesc)
    {
      if (current == NEWLINESYMBOL)
      {
        indesc = false;
      }
    } else
    {
      if (current == FASTASEPARATOR)
      {
        indesc = true;
      } else
      {
        if (!isspace((Ctypeargumenttype) current))
        {
          countcharacters++;
          switch (current)
          {
            case 'L':
            case 'I':
            case 'F':
            case 'E':
            case 'Q':
            case 'P':
            case 'X':
            case 'Z':
              countnonbases++;
              break;
            default:
              break;
          }
        }
      }
    }
  }
  if (countnonbases >= leastnumofnonbases)
  {
    return true;
  }
  return false;
}
*/

/* from guessprot.c

bool guessifproteinsequencestream2(const char *inputfile)
{
  Fgetcreturntype currentchar;
  Seqpos countnonbases = 0,
         countcharacters = 0;
  bool indesc = false;
  Genericstream inputstream;

  opengenericstream(&inputstream,inputfile);
  for (;;)
  {
    if (inputstream.isgzippedstream)
    {
      currentchar = gzgetc(inputstream.stream.gzippedstream);
    } else
    {
      currentchar = fgetc(inputstream.stream.fopenstream);
    }
    if (indesc)
    {
      if (currentchar == NEWLINESYMBOL)
      {
        indesc = false;
      }
    } else
    {
      if (currentchar == FASTASEPARATOR)
      {
        indesc = true;
      } else
      {
        if (!isspace((Ctypeargumenttype) currentchar))
        {
          countcharacters++;
          switch (currentchar)
          {
            case 'L':
            case 'I':
            case 'F':
            case 'E':
            case 'Q':
            case 'P':
            case 'X':
            case 'Z':
              countnonbases++;
              break;
            default:
              break;
          }
        }
      }
    }
    if (countcharacters >= (Seqpos) 1000)
    {
      break;
    }
  }
  closegenericstream(&inputstream,inputfile);
  if (countnonbases > 0 && countnonbases >= countcharacters/10)
  {
    return true;
  }
  return false;
}
*/
