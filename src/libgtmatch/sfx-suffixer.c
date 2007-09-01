/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
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
#define MAXPREFIXLENGTH ((((unsigned int) 1) << PREFIXLENBITS) - 1)
#define MAXCODEVALUE    ((((unsigned int) 1) << CODEBITS) - 1)

typedef struct
{
  unsigned int maxprefixlen:PREFIXLENBITS;
  unsigned int code:CODEBITS;
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
  unsigned int *filltable,
               *basepower;
  Seqpos *leftborder,
         *countspecialcodes,
         *suftab,
         *suftabptr;
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
  Suftabparts *suftabparts;
  const Encodedsequence *encseq;
  Readmode readmode;
  Seqpos widthofpart;
  unsigned int part,
               numofchars,
               prefixlength;
} Collectedsuffixes;

DECLAREARRAYSTRUCT(Seqpos);

static int initbasepower(unsigned int **basepower,
                         unsigned int **filltable,
                         unsigned int base,
                         unsigned int len,
                         Env *env)
{
  unsigned int thepower = (unsigned int) 1, i, minfailure;
  bool haserr = false;

  env_error_check(env);
  ALLOCASSIGNSPACE(*basepower,NULL,unsigned int,len+1);
  ALLOCASSIGNSPACE(*filltable,NULL,unsigned int,len);
  minfailure = UINT_MAX/base;
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
      if (firstspecial->specialpos > 0)
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

static Codetype codedownscale(const unsigned int *filltable,
                              const unsigned int *basepower,
                              Codetype code,
                              unsigned int prefixindex,
                              unsigned int maxprefixlen)
{
  unsigned int remain;

  code -= filltable[maxprefixlen];
  remain = maxprefixlen-prefixindex;
  code %= (filltable[remain]+1);
  code *= basepower[remain];
  code += filltable[prefixindex];
  return code;
}

static void derivespecialcodes(/*@unused@*/ const Encodedsequence *encseq,
                               Collectedsuffixes *csf,
                               bool deletevalues,
                               /*@unused@*/ Env *env)
{
  Codetype code;
  unsigned int prefixindex;
  unsigned long insertindex, j;
  Seqpos stidx;

  for (prefixindex=0; prefixindex < csf->prefixlength; prefixindex++)
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
          if (code != csf->filltable[0] || prefixindex > 0)
          {
            csf->countspecialcodes[FROMCODE2SPECIALCODE(code,
                                                        csf->numofchars)]++;
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
        if (prefixindex < csf->prefixlength - 1 &&
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

static int flushsuffixstarts(ArraySeqpos *suffixstarts,
                             Readmode readmode,
                             int(*processsuftab)(void *,const Seqpos *,
                                                 Readmode,Seqpos,Env *),
                             void *processsuftabinfo,
                             Env *env)
{
  if (suffixstarts->nextfreeSeqpos >= suffixstarts->allocatedSeqpos)
  {
    if (processsuftab(processsuftabinfo,
                      suffixstarts->spaceSeqpos,
                      readmode,
                      (Seqpos) suffixstarts->nextfreeSeqpos,
                      env) != 0)
    {
      return -1;
    }
    suffixstarts->nextfreeSeqpos = 0;
  }
  return 0;
}

static int insertfullspecialrange(ArraySeqpos *suffixstarts,
                                  Seqpos totallength,
                                  Readmode readmode,
                                  int(*processsuftab)(void *,const Seqpos *,
                                                      Readmode,Seqpos,Env *),
                                  void *processsuftabinfo,
                                  Seqpos leftpos,
                                  Seqpos rightpos,
                                  Env *env)
{
  Seqpos pos, tmppos;

  assert(leftpos < rightpos);
  /*
  printf("range %u,%u\n",leftpos,rightpos);
  */
  if(ISDIRREVERSE(readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while(true)
  {
    if(flushsuffixstarts(suffixstarts,readmode,processsuftab,
                         processsuftabinfo,env) != 0)
    {
      return -1;
    }
    if(ISDIRREVERSE(readmode))
    {
      tmppos = REVERSEPOS(totallength,pos);
      suffixstarts->spaceSeqpos[suffixstarts->nextfreeSeqpos++] 
        = tmppos;
      /*
      printf("map %u -> %u\n",pos,tmppos);
      */
      if(pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      suffixstarts->spaceSeqpos[suffixstarts->nextfreeSeqpos++] = pos;
      if(pos == rightpos-1)
      {
        break;
      }
      pos++;
    }
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
  ArraySeqpos suffixstarts;
  Sequencerange range;
  Seqpos totallength;
  bool haserr = false;

  totallength = getencseqtotallength(encseq);
  suffixstarts.spaceSeqpos = suftab;
  suffixstarts.allocatedSeqpos 
    = CALLCASTFUNC(Seqpos,unsigned_long,largestwidth);
  suffixstarts.nextfreeSeqpos = 0;
  if(hasspecialranges(encseq))
  {
    Specialrangeiterator *sri;

    sri = newspecialrangeiterator(encseq,ISDIRREVERSE(readmode) ? false : true,
                                  env);
    while(nextspecialrangeiterator(&range,sri))
    {
      if (insertfullspecialrange(&suffixstarts,
                                 totallength,
                                 readmode,
                                 processsuftab,
                                 processsuftabinfo,
                                 range.leftpos,
                                 range.rightpos,
                                 env) != 0)
      {
        haserr = true;
        break;
      }
    }
    freespecialrangeiterator(&sri,env);
  }
  if(!haserr)
  {
    if(flushsuffixstarts(&suffixstarts,readmode,processsuftab,
                         processsuftabinfo,env) != 0)
    {
      haserr = true;
    }
  }
  if(!haserr)
  {
    suffixstarts.spaceSeqpos[suffixstarts.nextfreeSeqpos++] = totallength;
    if (suffixstarts.nextfreeSeqpos > 0 && processsuftab != NULL)
    {
      if (processsuftab(processsuftabinfo,
                            suffixstarts.spaceSeqpos,
                            readmode,
                            (Seqpos) suffixstarts.nextfreeSeqpos,
                            env) != 0)
      {
        haserr = true;
      }
    }
  }
  return 0;
}

static int initsuffixerator(Collectedsuffixes *csf,
                            Seqpos specialcharacters,
                            Seqpos specialranges,
                            const Encodedsequence *encseq,
                            Readmode readmode,
                            unsigned int numofchars,
                            unsigned int prefixlength,
                            unsigned int numofparts,
                            Measuretime *mtime,
                            Env *env)
{
  unsigned int numofallcodes = 0, numofspecialcodes;
  Seqpos *optr;
  bool haserr = false;

  env_error_check(env);
  ALLOCASSIGNSPACE(csf->spaceCodeatposition,NULL,
                   Codeatposition,specialranges+1);
  csf->nextfreeCodeatposition = 0;
  csf->filltable = NULL;
  csf->basepower = NULL;
  csf->leftborder = NULL;
  csf->countspecialcodes = NULL;
  csf->suftab = NULL;
  csf->encseq = encseq;
  csf->readmode = readmode;
  csf->numofchars = numofchars;
  csf->prefixlength = prefixlength;
  csf->part = 0;
  if (prefixlength == 0 || prefixlength > MAXPREFIXLENGTH)
  {
    env_error_set(env,"argument for option -pl must be in the range [1,%u]",
                  MAXPREFIXLENGTH);
    haserr = true;
  }
  if (!haserr && initbasepower(&csf->basepower,
                               &csf->filltable,
                               numofchars,
                               prefixlength,
                               env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    assert(csf->basepower != NULL);
    numofallcodes = csf->basepower[prefixlength];
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
    Seqpos largestwidth;

    assert(csf->basepower != NULL);
    numofspecialcodes = csf->basepower[prefixlength-1];
    ALLOCASSIGNSPACE(csf->leftborder,NULL,Seqpos,numofallcodes+1);
    memset(csf->leftborder,0,
           sizeof (*csf->leftborder) * (size_t) numofallcodes);
    ALLOCASSIGNSPACE(csf->countspecialcodes,NULL,Seqpos,numofspecialcodes);
    memset(csf->countspecialcodes,0,
           sizeof (*csf->countspecialcodes) *
                  (size_t) numofspecialcodes);
    csf->storespecials = true;
    if (mtime != NULL)
    {
      deliverthetime(stdout,mtime,"counting prefix distribution",env);
    }
    getencseqkmers(encseq,
                   readmode,
                   updatekmercount,
                   csf,
                   numofchars,
                   prefixlength,
                   env);
    assert(specialranges+1 >= (Seqpos) csf->nextfreeCodeatposition);
    assert(csf->filltable != NULL);
    assert(csf->leftborder != NULL);
    /* printf("leftborder[0]=%u\n",csf.leftborder[0]); */
    for (optr = csf->leftborder + 1;
         optr < csf->leftborder + numofallcodes; optr++)
    {
      /*/ printf("leftborder[%u]=%u\n",(unsigned int) (optr - csf->leftborder),
                                   *optr); */
      *optr += *(optr-1);
    }
    csf->leftborder[numofallcodes]
      = getencseqtotallength(encseq) - specialcharacters;
    csf->suftabparts = newsuftabparts(numofparts,
                                      csf->leftborder,
                                      numofallcodes,
                                      getencseqtotallength(encseq)
                                        - specialcharacters,
                                      specialcharacters + 1,
                                      env);
    assert(csf->suftabparts != NULL);
    largestwidth = stpgetlargestwidth(csf->suftabparts);
    ALLOCASSIGNSPACE(csf->suftab,NULL,Seqpos,largestwidth);
    reversespecialcodes(csf->spaceCodeatposition,csf->nextfreeCodeatposition);
  }
  return haserr ? -1 : 0;
}

static bool preparethispart(Collectedsuffixes *csf,
                            Measuretime *mtime,
                            Env *env)
{
  if(csf->part >= stpgetnumofparts(csf->suftabparts))
  {
    return false;
  }
  csf->currentmincode = stpgetcurrentmincode(csf->part,csf->suftabparts);
  csf->currentmaxcode = stpgetcurrentmaxcode(csf->part,csf->suftabparts);
  csf->widthofpart = stpgetcurrentwidthofpart(csf->part,csf->suftabparts);
  csf->suftabptr = csf->suftab -
                   stpgetcurrentsuftaboffset(csf->part,csf->suftabparts);
  derivespecialcodes(NULL, /* not needed her */
                     csf,
                     (stpgetnumofparts(csf->suftabparts) == (unsigned int) 1)
                       ? true : false,
                     env);
  if (mtime != NULL)
  {
    deliverthetime(stdout,mtime,"inserting suffixes into buckets",env);
  }
  getencseqkmers(csf->encseq,
                 csf->readmode,
                 insertwithoutspecial,
                 csf,
                 csf->numofchars,
                 csf->prefixlength,
                 env);
  if (mtime != NULL)
  {
    deliverthetime(stdout,mtime,"sorting the buckets",env);
  }
  sortallbuckets(csf->suftabptr,
                 csf->encseq,
                 csf->readmode,
                 csf->leftborder,
                 csf->countspecialcodes,
                 csf->numofchars,
                 csf->prefixlength,
                 csf->currentmincode,
                 csf->currentmaxcode,
                 stpgetcurrentsumofwdith(csf->part,csf->suftabparts),
                 env);
  csf->part++;
  return true;
}

static void freeCollectedsuffixes(Collectedsuffixes *csf,Env *env)
{
  FREESPACE(csf->spaceCodeatposition);
  FREESPACE(csf->filltable);
  FREESPACE(csf->basepower);
  FREESPACE(csf->leftborder);
  FREESPACE(csf->countspecialcodes);
  FREESPACE(csf->suftab);
  freesuftabparts(csf->suftabparts,env);
}

int suffixerator(int(*processsuftab)(void *,const Seqpos *,
                                     Readmode,Seqpos,Env *),
                 void *processsuftabinfo,
                 Seqpos specialcharacters,
                 Seqpos specialranges,
                 const Encodedsequence *encseq,
                 Readmode readmode,
                 unsigned int numofchars,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 Measuretime *mtime,
                 Env *env)
{
  Collectedsuffixes csf;
  bool haserr = false;

  env_error_check(env);
  if(initsuffixerator(&csf,
                      specialcharacters,
                      specialranges,
                      encseq,
                      readmode,
                      numofchars,
                      prefixlength,
                      numofparts,
                      mtime,
                      env) != 0)
  {
    haserr = true;
  }
  if(!haserr)
  {
    Codetype specialcode;

    while(preparethispart(&csf,mtime,env))
    {
      if (processsuftab != NULL)
      {
        if (processsuftab(processsuftabinfo,
                          csf.suftab,
                          csf.readmode,
                          csf.widthofpart,
                          env) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
    if (insertallfullspecials(csf.encseq,
                              csf.readmode,
                              stpgetlargestwidth(csf.suftabparts),
                              csf.suftab,
                              processsuftab,
                              processsuftabinfo,
                              env) != 0)
    {
      haserr = true;
    }
    specialcode = FROMCODE2SPECIALCODE(csf.filltable[0],csf.numofchars);
    csf.countspecialcodes[specialcode] += (specialcharacters + 1);
  }
  freeCollectedsuffixes(&csf,env);
  return haserr ? -1 : 0;
}
