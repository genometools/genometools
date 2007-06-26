/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include "libgtcore/env.h"
#include "types.h"
#include "spacedef.h"
#include "arraydef.h"
#include "intbits.h"
#include "divmodmul.h"
#include "codespec.h"
#include "measure-time-if.h"
#include "partssuf-def.h"
#include "intcode-def.h"
#include "encseq-def.h"

#include "measure-time.pr"
#include "mappedstr.pr"
#include "partssuf.pr"
#include "bentleysedgewick.pr"

#define PREFIXLENBITS   4
#define CODEBITS        (32-PREFIXLENBITS)
#define MAXPREFIXLENGTH ((((int32_t) 1) << PREFIXLENBITS) - 1)
#define MAXCODEVALUE    ((((int32_t) 1) << CODEBITS) - 1)

typedef struct
{
  int32_t maxprefixlen:PREFIXLENBITS;
  int32_t code:CODEBITS;
  Seqpos position; /* get rid of this by using information from encseq */
} Codeatposition;

DECLAREARRAYSTRUCT(Codeatposition);

#undef LONGOUTPUT

typedef struct
{
  bool storespecials;
  Codetype currentmincode,
           currentmaxcode;
  int32_t *filltable,
          *basepower;
  Seqpos *leftborder,
         *countspecialcodes,
         *suftab,
         *suftabptr;
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
} Collectedsuffixes;

typedef struct
{
  Seqpos totallength,
         *spacesuffixstarts;
  unsigned long nextfreeUint, allocatedUint; 
  int(*processsuftab)(void *,const Seqpos *,Seqpos,Env *);
  void *processsuftabinfo;
  Env *env;
} InsertCompletespecials;


static int initbasepower(int32_t **basepower,
                         int32_t **filltable,
                         int32_t base,
                         int32_t len,
                         Env *env)
{
  int32_t thepower = (int32_t) 1, i, minfailure;
  bool haserr = false;

  env_error_check(env);
  ALLOCASSIGNSPACE(*basepower,NULL,int32_t,len+1);
  ALLOCASSIGNSPACE(*filltable,NULL,int32_t,len);
  minfailure = INT32_MAX/base;
  for (i=0; /* Nothing */; i++)
  {
    (*basepower)[i] = thepower;
    if (i == len)
    {
      break;
    }
    if (thepower >= minfailure)
    {
      env_error_set(env,"overflow when computing %u * %u\n",thepower,base);
      haserr = true;
      break;
    }
    thepower *= base;
  }
  if (!haserr)
  {
    for (i=0; i<len; i++)
    {
      (*filltable)[i] = (*basepower)[len-i]-1;
    }
  }
  if (haserr)
  {
    FREESPACE(*basepower);
    FREESPACE(*filltable);
    return -1;
  }
  return 0;
}

static void updatekmercount(void *processinfo,
                            Codetype code,
                            Seqpos position,
                            const Firstspecialpos *firstspecial,
                            Env *env)
{
  Collectedsuffixes *csf = (Collectedsuffixes *) processinfo;

  if (firstspecial->defined)
  {
    if (csf->storespecials)
    {
      env_error_check(env);

      if (firstspecial->specialpos > 0)
      {
        Codeatposition cp;

        cp.code = code;
        cp.maxprefixlen = firstspecial->specialpos;
        cp.position = (Seqpos) (position + firstspecial->specialpos);
        csf->spaceCodeatposition[csf->nextfreeCodeatposition++] = cp;
        csf->storespecials = false;
        csf->leftborder[code]++;
      }
    } else
    {
      if(firstspecial->specialpos > 0)
      {
        csf->leftborder[code]++;
      } else
      {
        csf->storespecials = true;
      }
    }
  } else
  {
    csf->leftborder[code]++;
  }
}

static void insertwithoutspecial(void *processinfo,
                                 Codetype code,
                                 Seqpos position,
                                 const Firstspecialpos *firstspecial,
                                 /*@unused@*/ Env *env)
{
  if (!firstspecial->defined)
  {
    Collectedsuffixes *csf = (Collectedsuffixes *) processinfo;

    if (code >= csf->currentmincode && code <= csf->currentmaxcode)
    {
      Seqpos stidx;

      stidx = --csf->leftborder[code];
#ifdef LONGOUTPUT
      printf("insert suffix %lu at location %lu\n",(Showuint) position,
                                                   (Showuint) stidx);
#endif
      csf->suftabptr[stidx] = position;
    }
  }
}

static void reversespecialcodes(Codeatposition *spaceCodeatposition,
                                unsigned long nextfreeCodeatposition)
{
  Codeatposition *front, *back, tmp;

  for (front = spaceCodeatposition,
       back = spaceCodeatposition + nextfreeCodeatposition - 1;
       front < back; front++, back--)
  {
    tmp = *front;
    *front = *back;
    *back = tmp;
  }
}

static Codetype codedownscale(const int32_t *filltable,
                              const int32_t *basepower,
                              Codetype code,
                              int32_t prefixindex,
                              int32_t maxprefixlen)
{
  int32_t remain;

  code -= filltable[maxprefixlen];
  remain = maxprefixlen-prefixindex;
  code %= (filltable[remain]+1);
  code *= basepower[remain];
  code += filltable[prefixindex];
  return code;
}

static void derivespecialcodes(/*@unused@*/ const Encodedsequence *encseq,
                               /*@unused@*/ Seqpos totallength,
                               unsigned int numofchars,
                               unsigned int prefixlength,
                               Collectedsuffixes *csf,
                               bool deletevalues,
                               /*@unused@*/ Env *env)
{
  Codetype code;
  unsigned int prefixindex;
  unsigned long insertindex, j, stidx;

  for (prefixindex=0; prefixindex < prefixlength; prefixindex++)
  {
    for (j=0, insertindex = 0; j < csf->nextfreeCodeatposition; j++)
    {
      if (prefixindex <= csf->spaceCodeatposition[j].maxprefixlen)
      {
        code = codedownscale(csf->filltable,
                             csf->basepower,
                             csf->spaceCodeatposition[j].code,
                             prefixindex,
                             csf->spaceCodeatposition[j].maxprefixlen);
        if (code >= csf->currentmincode &&
            code <= csf->currentmaxcode)
        {
          if(code != csf->filltable[0] || prefixindex > 0)
          {
            csf->countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)]++;
            stidx = --csf->leftborder[code];
#ifdef LONGOUTPUT
            printf("insert special_suffix %lu (code %u) at location %lu\n",
                   (Showuint) csf->spaceCodeatposition[j].position - prefixindex,
                   code,
                   (Showuint) stidx);
#endif
            csf->suftabptr[stidx] = csf->spaceCodeatposition[j].position -
                                    prefixindex;
          }
        }
      }
      if (deletevalues)
      {
        if (prefixindex < prefixlength - 1 &&
            prefixindex < csf->spaceCodeatposition[j].maxprefixlen)
        {
          if (insertindex < j)
          {
            csf->spaceCodeatposition[insertindex] =
              csf->spaceCodeatposition[j];
          }
          insertindex++;
        }
      }
    }
    if (deletevalues)
    {
      csf->nextfreeCodeatposition = insertindex;
    }
  }
}

static int insertfullspecialrange(InsertCompletespecials *ics,
                                  Seqpos startpos,
                                  Seqpos endpos)
{
  Seqpos pos;

  for (pos = startpos; pos < endpos; pos++)
  {
    if(ics->nextfreeUint >= ics->allocatedUint)
    {
      if (ics->processsuftab(ics->processsuftabinfo,
                             ics->spacesuffixstarts,
                             (Seqpos) ics->nextfreeUint,
                             ics->env) != 0)
      {
        return -1;
      }
      ics->nextfreeUint = 0;
    }
    ics->spacesuffixstarts[ics->nextfreeUint++] = pos;
  }
  return 0;
}

static int insertfullspecialpair(void *info,
                                 const PairSeqpos *pair,/*@unused@*/ Env *env)
{
  InsertCompletespecials *ics = (InsertCompletespecials *) info;

  if(insertfullspecialrange(ics,pair->uint0,pair->uint1) != 0)
  {
    return -1;
  }
  return 0;
}

static int insertallfullspecials(
                const Encodedsequence *encseq,
                Seqpos totallength,
                unsigned long largestwidth,
                Seqpos *suftab,
                int(*processsuftab)(void *,const Seqpos *,Seqpos,Env *),
                void *processsuftabinfo,
                Env *env)
{
  InsertCompletespecials ics;

  ics.totallength = totallength;
  ics.spacesuffixstarts = suftab;
  ics.allocatedUint = largestwidth;
  ics.nextfreeUint = 0;
  ics.processsuftab = processsuftab;
  ics.processsuftabinfo = processsuftabinfo;
  ics.env = env;
  if(overallspecialranges(encseq,insertfullspecialpair,&ics,env) != 0)
  {
    return -1;
  }
  if(insertfullspecialrange(&ics,totallength,totallength+1) != 0)
  {
    return -1;
  }
  if(ics.nextfreeUint > 0)
  {
    if (ics.processsuftab(ics.processsuftabinfo,
                          ics.spacesuffixstarts,
                          (Seqpos) ics.nextfreeUint,
                          env) != 0)
    {
      return -1;
    }
  }
  return 0;
}

int suffixerator(int(*processsuftab)(void *,const Seqpos *,Seqpos,Env *),
                 void *processsuftabinfo,
                 Seqpos totallength,
                 Seqpos specialcharacters,
                 Seqpos specialranges,
                 const Encodedsequence *encseq,
                 unsigned int numofchars,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 Measuretime *mtime,
                 Env *env)
{
  unsigned int numofallcodes = 0, numofspecialcodes, part;
  Seqpos *optr;
  Collectedsuffixes csf;
  Suftabparts *suftabparts = NULL;
  bool haserr = false;

  env_error_check(env);
  ALLOCASSIGNSPACE(csf.spaceCodeatposition,NULL,
                   Codeatposition,specialranges+1);
  csf.nextfreeCodeatposition = 0;
  csf.filltable = NULL;
  csf.basepower = NULL;
  csf.leftborder = NULL;
  csf.countspecialcodes = NULL;
  csf.suftab = NULL;
  if (prefixlength == 0 || prefixlength > MAXPREFIXLENGTH)
  {
    env_error_set(env,"argument for option -pl must be in the range [1,%u]",
                  MAXPREFIXLENGTH);
    haserr = true;
  }
  if (!haserr && initbasepower(&csf.basepower,
                               &csf.filltable,
                               numofchars,
                               prefixlength,
                               env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    assert(csf.basepower != NULL);
    numofallcodes = csf.basepower[prefixlength];
    if (numofallcodes-1 > MAXCODEVALUE)
    {
      env_error_set(env,
                    "alphasize^prefixlength-1 = %u does not fit into "
                    " %lu bits: choose smaller value for prefixlength",
                    numofallcodes-1,
                    (Showuint) CODEBITS);
      haserr = true;
    }
  }
  if (!haserr)
  {
    Codetype specialcode;
    Seqpos largestwidth;

    assert(csf.basepower != NULL);
    numofspecialcodes = csf.basepower[prefixlength-1];
    ALLOCASSIGNSPACE(csf.leftborder,NULL,Seqpos,numofallcodes+1);
    memset(csf.leftborder,0,
           sizeof (*csf.leftborder) * (size_t) numofallcodes);
    ALLOCASSIGNSPACE(csf.countspecialcodes,NULL,Seqpos,numofspecialcodes);
    memset(csf.countspecialcodes,0,
           sizeof (*csf.countspecialcodes) *
                  (size_t) numofspecialcodes);
    csf.storespecials = true;
    deliverthetime(stdout,mtime,"counting prefix distribution",env);
    getencseqkmers(encseq,
                   totallength,
                   updatekmercount,
                   &csf,
                   numofchars,
                   prefixlength,
                   env);
    assert(csf.nextfreeCodeatposition <= specialranges+1);
    assert(csf.filltable != NULL);
    assert(csf.leftborder != NULL);
    for (optr = csf.leftborder + 1;
         optr < csf.leftborder + numofallcodes; optr++)
    {
      *optr += *(optr-1);
    }
    csf.leftborder[numofallcodes] = (Seqpos) (totallength-specialcharacters);
    suftabparts = newsuftabparts(numofparts,
                                 csf.leftborder,
                                 numofallcodes,
                                 totallength - specialcharacters,
                                 specialcharacters + 1,
                                 env);
    assert(suftabparts != NULL);
    largestwidth = stpgetlargestwidth(suftabparts);
    ALLOCASSIGNSPACE(csf.suftab,NULL,Seqpos,largestwidth);
    reversespecialcodes(csf.spaceCodeatposition,csf.nextfreeCodeatposition);
    for (part = 0; part < stpgetnumofparts(suftabparts); part++)
    {
      csf.currentmincode = stpgetcurrentmincode(part,suftabparts);
      csf.suftabptr = csf.suftab - stpgetcurrentsuftaboffset(part,suftabparts);
      csf.currentmaxcode = stpgetcurrentmaxcode(part,suftabparts);
      derivespecialcodes(encseq,
                         totallength,
                         numofchars,
                         prefixlength,
                         &csf,
                         (stpgetnumofparts(suftabparts) == (unsigned int) 1)
                           ? true : false,
                         env);
      deliverthetime(stdout,mtime,"inserting suffixes into buckets",env);
      getencseqkmers(encseq,
                     totallength,
                     insertwithoutspecial,
                     &csf,
                     numofchars,
                     prefixlength,
                     env);
      deliverthetime(stdout,mtime,"sorting the buckets",env);
      sortallbuckets(csf.suftabptr,
                     encseq,
                     csf.leftborder,
                     csf.countspecialcodes,
                     totallength,
                     numofchars,
                     prefixlength,
                     csf.currentmincode,
                     csf.currentmaxcode,
                     stpgetcurrentsumofwdith(part,suftabparts),
                     env);
      if (processsuftab != NULL)
      {
        if (processsuftab(processsuftabinfo,
                          csf.suftab,
                          stpgetcurrentwidtofpart(part,suftabparts),
                          env) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
    if (insertallfullspecials(encseq,
                              totallength,
                              largestwidth,
                              csf.suftab,
                              processsuftab,
                              processsuftabinfo,
                              env) != 0)
    {
      haserr = true;
    }
    specialcode = FROMCODE2SPECIALCODE(csf.filltable[0],numofchars);
    csf.countspecialcodes[specialcode] += (specialcharacters + 1);
  }
  FREESPACE(csf.spaceCodeatposition);
  FREESPACE(csf.filltable);
  FREESPACE(csf.basepower);
  FREESPACE(csf.leftborder);
  FREESPACE(csf.countspecialcodes);
  FREESPACE(csf.suftab);
  freesuftabparts(suftabparts,env);
  return haserr ? -1 : 0;
}
