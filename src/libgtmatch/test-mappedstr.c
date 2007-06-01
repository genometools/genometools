/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "types.h"
#include "arraydef.h"
#include "spacedef.h"
#include "alphadef.h"
#include "chardef.h"
#include "addnextchar.h"
#include "encseq-def.h"
#include "sarr-def.h"

#include "kmer2string.pr"
#include "compfilenm.pr"
#include "mappedstr.pr"
#include "alphabet.pr"
#include "sfxmap.pr"

static Uint qgram2codefillspecial(Uint numofchars,
                                  unsigned int kmersize,
                                  const Encodedsequence *encseq,
                                  Uint64 startpos,
                                  Uint64 totallength)
{
  Uint integercode;
  Uint64 pos;
  bool foundspecial;
  Uchar cc;

  if (startpos >= totallength)
  {
    integercode = numofchars - 1;
    foundspecial = true;
  } else
  {
    cc = getencodedchar64(encseq,startpos);
    if (ISSPECIAL(cc))
    {
      integercode = numofchars - 1;
      foundspecial = true;
    } else
    {
      integercode = (Uint) cc;
      foundspecial = false;
    }
  }
  for (pos = startpos + (Uint64) 1; pos< startpos + (Uint64) kmersize; pos++)
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
        cc = getencodedchar64(encseq,pos);
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

static void outkmeroccurrence(void *processinfo,
                              Uint code,
                              /*@unused@*/ Uint64 position,
                              /*@unused@*/ const DefinedUint
                                                   *firstspecialposition,
                              Env *env)
{
  ArrayUint *codelist = (ArrayUint *) processinfo;

  STOREINARRAY(codelist,Uint,1024,code);
}

static void collectkmercode(ArrayUint *codelist,
                            const Encodedsequence *encseq,
                            unsigned int kmersize,
                            Uint numofchars,
                            Uint64 stringtotallength,
                            Env *env)
{
  Uint64 offset;
  Uint code;

  for (offset=0; offset<=stringtotallength; offset++)
  {
    code = qgram2codefillspecial(numofchars,
                                 kmersize,
                                 encseq,
                                 offset,
                                 stringtotallength);
    STOREINARRAY(codelist,Uint,1024,code);
  }
}

static int comparecodelists(const ArrayUint *codeliststream,
                            const ArrayUint *codeliststring,
                            unsigned int kmersize,
                            Uint numofchars,
                            const char *characters,
                            Env *env)
{
  Uint i;
  char buffer1[64+1], buffer2[64+1];

  if (codeliststream->nextfreeUint != codeliststring->nextfreeUint)
  {
    env_error_set(env,
                  "length codeliststream= %lu != %lu =length codeliststring",
                  (Showuint) codeliststream->nextfreeUint,
                  (Showuint) codeliststring->nextfreeUint);
    return -1;
  }
  for (i=0; i<codeliststream->nextfreeUint; i++)
  {
    if (codeliststream->spaceUint[i] != codeliststring->spaceUint[i])
    {
      kmercode2string(buffer1,
                      codeliststream->spaceUint[i],
                      numofchars,
                      kmersize,
                      characters);
      kmercode2string(buffer2,
                      codeliststring->spaceUint[i],
                      numofchars,
                      kmersize,
                      characters);
      env_error_set(env,
                    "codeliststream[%lu] = %lu != %lu = "
                    "codeliststring[%lu]\n%s != %s",
                    (Showuint) i,
                    (Showuint) codeliststream->spaceUint[i],
                    (Showuint) codeliststring->spaceUint[i],
                    (Showuint) i,
                    buffer1,
                    buffer2);
      return -1;
    }
  }
  return 0;
}

static int verifycodelists(const Encodedsequence *encseq,
                           const Uchar *characters,
                           unsigned int kmersize,
                           Uint numofchars,
                           Uint64 stringtotallength,
                           const ArrayUint *codeliststream,
                           Env *env)
{
  bool haserr = false;
  ArrayUint codeliststring;

  INITARRAY(&codeliststring,Uint);
  collectkmercode(&codeliststring,
                  encseq,
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
  FREEARRAY(&codeliststring,Uint);
  return haserr ? -1 : 0;
}

int verifymappedstr(const Suffixarray *suffixarray,Env *env)
{
  Uint numofchars;
  ArrayUint codeliststream;
  bool haserr = false;

  numofchars = getnumofcharsAlphabet(suffixarray->alpha);
  INITARRAY(&codeliststream,Uint);
  if (getfastastreamkmers((const char **) suffixarray->filenametab,
                          suffixarray->numoffiles,
                          outkmeroccurrence,
                          &codeliststream,
                          numofchars,
                          suffixarray->prefixlength,
                          getsymbolmapAlphabet(suffixarray->alpha),
                          env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (verifycodelists(suffixarray->encseq,
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
  FREEARRAY(&codeliststream,Uint);
  return haserr ? -1 : 0;
}
