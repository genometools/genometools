/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "libgtcore/env.h"
#include "spacedef.h"
#include "arraydef.h"
#include "intbits.h"
#include "divmodmul.h"
#include "measure-time-if.h"
#include "intcode-def.h"
#include "encseq-def.h"
#include "safecast-gen.h"
#include "sfx-codespec.h"
#include "sfx-partssuf-def.h"

#include "measure-time.pr"
#include "sfx-mappedstr.pr"
#include "sfx-partssuf.pr"
#include "sfx-bentsedg.pr"

#define CODEBITS        (32-PREFIXLENBITS)
#define MAXPREFIXLENGTH ((((uint32_t) 1) << PREFIXLENBITS) - 1)
#define MAXCODEVALUE    ((((uint32_t) 1) << CODEBITS) - 1)

typedef struct
{
  uint32_t maxprefixlen:PREFIXLENBITS;
  uint32_t code:CODEBITS;
  Seqpos position; /* get rid of this by using information from encseq */
} Codeatposition;

DECLAREARRAYSTRUCT(Codeatposition);

#define LONGOUTPUT
#undef LONGOUTPUT

typedef struct
{
  bool storespecials;
  Codetype currentmincode,
           currentmaxcode;
  uint32_t *filltable,
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
  unsigned long nextfreeUint, 
                allocatedUint; 
  int(*processsuftab)(void *,const Seqpos *,Readmode,Seqpos,Env *);
  void *processsuftabinfo;
  Readmode readmode;
  Env *env;
} InsertCompletespecials;

static int initbasepower(uint32_t **basepower,
                         uint32_t **filltable,
                         uint32_t base,
                         uint32_t len,
                         Env *env)
{
  uint32_t thepower = (uint32_t) 1, i, minfailure;
  bool haserr = false;

  env_error_check(env);
  ALLOCASSIGNSPACE(*basepower,NULL,uint32_t,len+1);
  ALLOCASSIGNSPACE(*filltable,NULL,uint32_t,len);
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
      env_error_set(env,"overflow when computing %u * %u",thepower,base);
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
        cp.position = position + firstspecial->specialpos;
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
      printf("insert suffix " FormatSeqpos " at location " FormatSeqpos "\n",
              PRINTSeqposcast(position),
              PRINTSeqposcast(stidx));
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

static Codetype codedownscale(const uint32_t *filltable,
                              const uint32_t *basepower,
                              Codetype code,
                              uint32_t prefixindex,
                              uint32_t maxprefixlen)
{
  uint32_t remain;

  code -= filltable[maxprefixlen];
  remain = maxprefixlen-prefixindex;
  code %= (filltable[remain]+1);
  code *= basepower[remain];
  code += filltable[prefixindex];
  return code;
}

static void derivespecialcodes(/*@unused@*/ const Encodedsequence *encseq,
                               uint32_t numofchars,
                               uint32_t prefixlength,
                               Collectedsuffixes *csf,
                               bool deletevalues,
                               /*@unused@*/ Env *env)
{
  Codetype code;
  uint32_t prefixindex;
  unsigned long insertindex, j;
  Seqpos stidx;

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
            printf("insert special_suffix " FormatSeqpos 
                   " (code %u) at location " FormatSeqpos "\n",
                   PRINTSeqposcast(csf->spaceCodeatposition[j].position - 
                                   prefixindex),
                   (unsigned int) code,
                   PRINTSeqposcast(stidx));
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
                             ics->readmode,
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
                                 const Sequencerange *pair,
                                 /*@unused@*/ Env *env)
{
  InsertCompletespecials *ics = (InsertCompletespecials *) info;

  if(insertfullspecialrange(ics,pair->leftpos,pair->rightpos) != 0)
  {
    return -1;
  }
  return 0;
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

static int insertallfullspecials(
                const Encodedsequence *encseq,
                Readmode readmode,
                Seqpos largestwidth,
                Seqpos *suftab,
                int(*processsuftab)(void *,const Seqpos *,Readmode,
                                    Seqpos,Env *),
                void *processsuftabinfo,
                Env *env)
{
  InsertCompletespecials ics;

  ics.totallength = getencseqtotallength(encseq);
  ics.spacesuffixstarts = suftab;
  ics.allocatedUint = CALLCASTFUNC(Seqpos,unsigned_long,largestwidth);
  ics.nextfreeUint = 0;
  ics.processsuftab = processsuftab;
  ics.processsuftabinfo = processsuftabinfo;
  ics.readmode = readmode;
  ics.env = env;
  if(overallspecialranges(encseq,
                          readmode,
                          insertfullspecialpair,&ics,env) != 0)
  {
    return -1;
  }
  if(insertfullspecialrange(&ics,ics.totallength,ics.totallength+1) != 0)
  {
    return -2;
  }
  if(ics.nextfreeUint > 0)
  {
    if (ics.processsuftab(ics.processsuftabinfo,
                          ics.spacesuffixstarts,
                          ics.readmode,
                          (Seqpos) ics.nextfreeUint,
                          env) != 0)
    {
      return -3;
    }
  }
  return 0;
}

int suffixerator(int(*processsuftab)(void *,const Seqpos *,
                                     Readmode,Seqpos,Env *),
                 void *processsuftabinfo,
                 Seqpos specialcharacters,
                 Seqpos specialranges,
                 const Encodedsequence *encseq,
                 Readmode readmode,
                 uint32_t numofchars,
                 uint32_t prefixlength,
                 uint32_t numofparts,
                 Measuretime *mtime,
                 Env *env)
{
  uint32_t numofallcodes = 0, numofspecialcodes, part;
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
                    " %u bits: choose smaller value for prefixlength",
                    numofallcodes-1,
                    (unsigned int) CODEBITS);
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
                   readmode,
                   updatekmercount,
                   &csf,
                   numofchars,
                   prefixlength,
                   env);
    assert(specialranges+1 >= (Seqpos) csf.nextfreeCodeatposition);
    assert(csf.filltable != NULL);
    assert(csf.leftborder != NULL);
    // printf("leftborder[0]=%u\n",csf.leftborder[0]);
    for (optr = csf.leftborder + 1;
         optr < csf.leftborder + numofallcodes; optr++)
    {
      // printf("leftborder[%u]=%u\n",(unsigned int) (optr - csf.leftborder),
                                   // *optr);
      *optr += *(optr-1);
    }
    csf.leftborder[numofallcodes] 
      = getencseqtotallength(encseq) - specialcharacters;
    suftabparts = newsuftabparts(numofparts,
                                 csf.leftborder,
                                 numofallcodes,
                                 getencseqtotallength(encseq)
                                   - specialcharacters,
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
      derivespecialcodes(NULL, /* not needed her */
                         numofchars,
                         prefixlength,
                         &csf,
                         (stpgetnumofparts(suftabparts) == (uint32_t) 1)
                           ? true : false,
                         env);
      deliverthetime(stdout,mtime,"inserting suffixes into buckets",env);
      getencseqkmers(encseq,
                     readmode,
                     insertwithoutspecial,
                     &csf,
                     numofchars,
                     prefixlength,
                     env);
      deliverthetime(stdout,mtime,"sorting the buckets",env);
      sortallbuckets(csf.suftabptr,
                     encseq,
                     readmode,
                     csf.leftborder,
                     csf.countspecialcodes,
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
                          readmode,
                          stpgetcurrentwidtofpart(part,suftabparts),
                          env) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
    if (insertallfullspecials(encseq,
                              readmode,
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
