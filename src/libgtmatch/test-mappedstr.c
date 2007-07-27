/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "arraydef.h"
#include "spacedef.h"
#include "alphadef.h"
#include "chardef.h"
#include "encseq-def.h"
#include "intcode-def.h"
#include "sarr-def.h"
#include "sfx-nextchar.h"

#include "alphabet.pr"
#include "kmer2string.pr"
#include "sfx-mappedstr.pr"
#include "sfx-map.pr"

static Codetype qgram2codefillspecial(uint32_t numofchars,
                                      uint32_t kmersize,
                                      const Encodedsequence *encseq,
                                      Readmode readmode,
                                      Seqpos startpos,
                                      Seqpos totallength)
{
  Codetype integercode;
  Seqpos pos;
  bool foundspecial;
  Uchar cc;

  if (startpos >= totallength)
  {
    integercode = numofchars - 1;
    foundspecial = true;
  } else
  {
    cc = getencodedchar(encseq,startpos,readmode);
    if (ISSPECIAL(cc))
    {
      integercode = numofchars - 1;
      foundspecial = true;
    } else
    {
      integercode = (Codetype) cc;
      foundspecial = false;
    }
  }
  for (pos = startpos + 1; pos < startpos + kmersize; pos++)
  {
    if (foundspecial)
    {
      ADDNEXTCHAR(integercode,numofchars-1,numofchars);
    } else
    {
      if (pos >= totallength)
      {
        ADDNEXTCHAR(integercode,numofchars-1,numofchars);
        foundspecial = true;
      } else
      {
        cc = getencodedchar(encseq,pos,readmode);
        if (ISSPECIAL(cc))
        {
          ADDNEXTCHAR(integercode,numofchars-1,numofchars);
          foundspecial = true;
        } else
        {
          ADDNEXTCHAR(integercode,cc,numofchars);
        }
      }
    }
  }
  return integercode;
}

DECLAREARRAYSTRUCT(Codetype);

static void outkmeroccurrence(void *processinfo,
                              Codetype code,
                              /*@unused@*/ Seqpos position,
                              /*@unused@*/ const Firstspecialpos
                                                 *firstspecialposition,
                              Env *env)
{
  ArrayCodetype *codelist = (ArrayCodetype *) processinfo;

  STOREINARRAY(codelist,Codetype,1024,code);
}

static void collectkmercode(ArrayCodetype *codelist,
                            const Encodedsequence *encseq,
                            Readmode readmode,
                            uint32_t kmersize,
                            uint32_t numofchars,
                            Seqpos stringtotallength,
                            Env *env)
{
  Seqpos offset;
  Codetype code;

  for (offset=0; offset<=stringtotallength; offset++)
  {
    code = qgram2codefillspecial(numofchars,
                                 kmersize,
                                 encseq,
                                 readmode,
                                 offset,
                                 stringtotallength);
    STOREINARRAY(codelist,Codetype,1024,code);
  }
}

static int comparecodelists(const ArrayCodetype *codeliststream,
                            const ArrayCodetype *codeliststring,
                            uint32_t kmersize,
                            uint32_t numofchars,
                            const char *characters,
                            Env *env)
{
  unsigned long i;
  char buffer1[64+1], buffer2[64+1];

  if (codeliststream->nextfreeCodetype != codeliststring->nextfreeCodetype)
  {
    env_error_set(env,
                  "length codeliststream= %lu != %lu =length codeliststring",
                  (unsigned long) codeliststream->nextfreeCodetype,
                  (unsigned long) codeliststring->nextfreeCodetype);
    return -1;
  }
  for (i=0; i<codeliststream->nextfreeCodetype; i++)
  {
    if (codeliststream->spaceCodetype[i] != codeliststring->spaceCodetype[i])
    {
      kmercode2string(buffer1,
                      codeliststream->spaceCodetype[i],
                      numofchars,
                      kmersize,
                      characters);
      kmercode2string(buffer2,
                      codeliststring->spaceCodetype[i],
                      numofchars,
                      kmersize,
                      characters);
      env_error_set(env,
                    "codeliststream[%lu] = %u != %u = "
                    "codeliststring[%lu]\n%s != %s",
                    i,
                    codeliststream->spaceCodetype[i],
                    codeliststring->spaceCodetype[i],
                    i,
                    buffer1,
                    buffer2);
      return -1;
    }
  }
  return 0;
}

static int verifycodelists(const Encodedsequence *encseq,
                           Readmode readmode,
                           const Uchar *characters,
                           uint32_t kmersize,
                           uint32_t numofchars,
                           Seqpos stringtotallength,
                           const ArrayCodetype *codeliststream,
                           Env *env)
{
  bool haserr = false;
  ArrayCodetype codeliststring;

  INITARRAY(&codeliststring,Codetype);
  collectkmercode(&codeliststring,
                  encseq,
                  readmode,
                  kmersize,
                  numofchars,
                  stringtotallength,
                  env);
  if (comparecodelists(codeliststream,
                       &codeliststring,
                       kmersize,
                       numofchars,
                       (const char *) characters,
                       env) != 0)
  {
    haserr = true;
  }
  FREEARRAY(&codeliststring,Codetype);
  return haserr ? -1 : 0;
}

int verifymappedstr(const Suffixarray *suffixarray,Env *env)
{
  uint32_t numofchars;
  ArrayCodetype codeliststream;
  bool haserr = false;

  numofchars = getnumofcharsAlphabet(suffixarray->alpha);
  INITARRAY(&codeliststream,Codetype);
  if (getfastastreamkmers(suffixarray->filenametab,
                          outkmeroccurrence,
                          &codeliststream,
                          numofchars,
                          suffixarray->prefixlength,
                          getsymbolmapAlphabet(suffixarray->alpha),
                          false,
                          env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (verifycodelists(suffixarray->encseq,
                        suffixarray->readmode,
                        getcharactersAlphabet(suffixarray->alpha),
                        suffixarray->prefixlength,
                        numofchars,
                        getencseqtotallength(suffixarray->encseq),
                        &codeliststream,
                        env) != 0)
    {
      haserr = true;
    }
  }
  FREEARRAY(&codeliststream,Codetype);
  return haserr ? -1 : 0;
}
