/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
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

DECLAREARRAYSTRUCT(Seqpos);

typedef enum
{
  Partphase,
  Specialrangephase,
  Finalsuffixphase
} Suffixeratorphase;

typedef struct
{
  bool storespecials;
  Codetype currentmincode,
           currentmaxcode;
  unsigned int *filltable,
               *basepower;
  Seqpos specialcharacters,
         *leftborder,
         *countspecialcodes,
         *suftab,
         *suftabptr;
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
  Suftabparts *suftabparts;
  const Encodedsequence *encseq;
  Readmode readmode;
  Seqpos widthofpart;
  Seqpos totallength;
  unsigned int part,
               numofchars,
               prefixlength;
  ArraySeqpos fusp;
  Specialrangeiterator *sri;
  Suffixeratorphase phase;
  Sequencerange overhang;
  bool exhausted;
} Sfxiterator;

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
  Sfxiterator *csf = (Sfxiterator *) processinfo;

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
    Sfxiterator *csf = (Sfxiterator *) processinfo;

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
                               Sfxiterator *csf,
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

static int processfullspecialrange(Sfxiterator *csf,
                                   int(*processsuftab)(void *,const Seqpos *,
                                                       Readmode,Seqpos,Env *),
                                   void *processsuftabinfo,
                                   Seqpos leftpos,
                                   Seqpos rightpos,
                                   Env *env)
{
  Seqpos pos;
  bool stop = false;

  assert(leftpos < rightpos);
  if (ISDIRREVERSE(csf->readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while (!stop)
  {
    if (ISDIRREVERSE(csf->readmode))
    {
      csf->fusp.spaceSeqpos[csf->fusp.nextfreeSeqpos++]
        = REVERSEPOS(csf->totallength,pos);
      if (pos == leftpos)
      {
        stop = true;
      } else
      {
        pos--;
      }
    } else
    {
      csf->fusp.spaceSeqpos[csf->fusp.nextfreeSeqpos++] = pos;
      if (pos == rightpos-1)
      {
        stop = true;
      } else
      {
        pos++;
      }
    }
    if (csf->fusp.nextfreeSeqpos >= csf->fusp.allocatedSeqpos)
    {
      if (processsuftab(processsuftabinfo,
                        csf->fusp.spaceSeqpos,
                        csf->readmode,
                        (Seqpos) csf->fusp.nextfreeSeqpos,
                        env) != 0)
      {
        return -1;
      }
      csf->fusp.nextfreeSeqpos = 0;
    }
  }
  return 0;
}

static void insertfullspecialrange(Sfxiterator *csf,
                                   Seqpos leftpos,
                                   Seqpos rightpos)
{
  Seqpos pos;

  assert(leftpos < rightpos);
  if (ISDIRREVERSE(csf->readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while (true)
  {
    if (ISDIRREVERSE(csf->readmode))
    {
      csf->fusp.spaceSeqpos[csf->fusp.nextfreeSeqpos++]
        = REVERSEPOS(csf->totallength,pos);
      if (pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      csf->fusp.spaceSeqpos[csf->fusp.nextfreeSeqpos++] = pos;
      if (pos == rightpos-1)
      {
        break;
      }
      pos++;
    }
  }
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

static int initsuffixerator(Sfxiterator *csf,
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
  csf->totallength = getencseqtotallength(encseq);
  csf->specialcharacters = specialcharacters;
  csf->part = 0;
  csf->phase = Partphase;
  csf->exhausted = false;
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
    csf->leftborder[numofallcodes] = csf->totallength - specialcharacters;
    csf->suftabparts = newsuftabparts(numofparts,
                                      csf->leftborder,
                                      numofallcodes,
                                      csf->totallength - specialcharacters,
                                      specialcharacters + 1,
                                      env);
    assert(csf->suftabparts != NULL);
    ALLOCASSIGNSPACE(csf->suftab,NULL,Seqpos,
                     stpgetlargestwidth(csf->suftabparts));
    reversespecialcodes(csf->spaceCodeatposition,csf->nextfreeCodeatposition);
    /* XXX the following can be removed for the iterator version */
    csf->fusp.spaceSeqpos = csf->suftab;
    csf->fusp.allocatedSeqpos
    = CALLCASTFUNC(Seqpos,unsigned_long,
                   stpgetlargestwidth(csf->suftabparts));
    csf->fusp.nextfreeSeqpos = 0;
  }
  return haserr ? -1 : 0;
}

static bool preparethispart(Sfxiterator *csf,
                            Measuretime *mtime,
                            Env *env)
{
  if (csf->part >= stpgetnumofparts(csf->suftabparts))
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

static void freeSfxiterator(Sfxiterator *csf,Env *env)
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
  Sfxiterator csf;
  bool haserr = false;

  env_error_check(env);
  if (initsuffixerator(&csf,
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
  if (!haserr)
  {
    Codetype specialcode;

    while (preparethispart(&csf,mtime,env))  /* parts phase */
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
    if (!haserr && hasspecialranges(csf.encseq)) /* special ranges phase */
    {
      Sequencerange range;
      Specialrangeiterator *sri;

      sri = newspecialrangeiterator(csf.encseq,
                                    ISDIRREVERSE(csf.readmode) ? false : true,
                                    env);
      while (nextspecialrangeiterator(&range,sri))
      {
        if (processfullspecialrange(&csf,
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
    if (!haserr) /* last special phase */
    {
      csf.fusp.spaceSeqpos[csf.fusp.nextfreeSeqpos++] = csf.totallength;
      if (processsuftab != NULL)
      {
        if (processsuftab(processsuftabinfo,
                          csf.fusp.spaceSeqpos,
                          csf.readmode,
                          (Seqpos) csf.fusp.nextfreeSeqpos,
                          env) != 0)
        {
          haserr = true;
        }
      }
    }
    specialcode = FROMCODE2SPECIALCODE(csf.filltable[0],csf.numofchars);
    csf.countspecialcodes[specialcode] += (specialcharacters + 1);
  }
  freeSfxiterator(&csf,env);
  return haserr ? -1 : 0;
}

 void freeSfxiterator2(Sfxiterator **sfxiterator,Env *env)
{
  Codetype specialcode;

  specialcode = FROMCODE2SPECIALCODE((*sfxiterator)->filltable[0],
                                    (*sfxiterator)->numofchars);
  (*sfxiterator)->countspecialcodes[specialcode]
    += ((*sfxiterator)->specialcharacters + 1);
  freeSfxiterator(*sfxiterator,env);
  FREESPACE(*sfxiterator);
}

 Sfxiterator *newsfxiterator(Seqpos specialcharacters,
                             Seqpos specialranges,
                             const Encodedsequence *encseq,
                             Readmode readmode,
                             unsigned int numofchars,
                             unsigned int prefixlength,
                             unsigned int numofparts,
                             Measuretime *mtime,
                             Env *env)
{
  Sfxiterator *csf;

  env_error_check(env);
  ALLOCASSIGNSPACE(csf,NULL,Sfxiterator,1);
  if (initsuffixerator(csf,
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
    freeSfxiterator2(&csf,env);
    return NULL;
  }
  return csf;
}

static void fillspecialnextpage(Sfxiterator *csf)
{
  Sequencerange range;
  Seqpos width;

  while (true)
  {
    if (csf->overhang.leftpos < csf->overhang.rightpos)
    {
      width = csf->overhang.rightpos - csf->overhang.leftpos;
      if (csf->fusp.nextfreeSeqpos + width > csf->fusp.allocatedSeqpos)
      {
        /* does not fit into the buffer, so only output a part */
        unsigned long rest = csf->fusp.nextfreeSeqpos +
                             width - csf->fusp.allocatedSeqpos;
        assert(rest > 0);
        if (ISDIRREVERSE(csf->readmode))
        {
          insertfullspecialrange(csf,csf->overhang.leftpos + rest,
                                 csf->overhang.rightpos);
          csf->overhang.rightpos = csf->overhang.leftpos + rest;
        } else
        {
          insertfullspecialrange(csf,csf->overhang.leftpos,
                                     csf->overhang.rightpos - rest);
          csf->overhang.leftpos = csf->overhang.rightpos - rest;
        }
        break;
      }
      if (csf->fusp.nextfreeSeqpos + width == csf->fusp.allocatedSeqpos)
      { /* overhang fits into the buffer and buffer is full */
        insertfullspecialrange(csf,csf->overhang.leftpos,
                               csf->overhang.rightpos);
        csf->overhang.leftpos = csf->overhang.rightpos = 0;
        break;
      }
      /* overhang fits into the buffer and buffer is not full */
      insertfullspecialrange(csf,csf->overhang.leftpos,
                             csf->overhang.rightpos);
      csf->overhang.leftpos = csf->overhang.rightpos = 0;
    } else
    {
      if (nextspecialrangeiterator(&range,csf->sri))
      {
        width = range.rightpos - range.leftpos;
        if (csf->fusp.nextfreeSeqpos + width > csf->fusp.allocatedSeqpos)
        { /* does not fit into the buffer, so only output a part */
          unsigned long rest = csf->fusp.nextfreeSeqpos +
                               width - csf->fusp.allocatedSeqpos;
          if (ISDIRREVERSE(csf->readmode))
          {
            insertfullspecialrange(csf,range.leftpos + rest,
                                   range.rightpos);
            csf->overhang.leftpos = range.leftpos;
            csf->overhang.rightpos = range.leftpos + rest;
          } else
          {
            insertfullspecialrange(csf,range.leftpos,range.rightpos - rest);
            csf->overhang.leftpos = range.rightpos - rest;
            csf->overhang.rightpos = range.rightpos;
          }
          break;
        }
        if (csf->fusp.nextfreeSeqpos + width == csf->fusp.allocatedSeqpos)
        { /* overhang fits into the buffer and buffer is full */
          insertfullspecialrange(csf,range.leftpos,range.rightpos);
          csf->overhang.leftpos = csf->overhang.rightpos = 0;
          break;
        }
        insertfullspecialrange(csf,range.leftpos,range.rightpos);
        csf->overhang.leftpos = csf->overhang.rightpos = 0;
      } else
      {
        if(csf->fusp.nextfreeSeqpos < csf->fusp.allocatedSeqpos)
        {
          csf->fusp.spaceSeqpos[csf->fusp.nextfreeSeqpos++] = csf->totallength;
          csf->exhausted = true;
        }
        break;
      }
    }
  }
}

 const Seqpos *nextSfxiterator(Seqpos *len,Measuretime *mtime,
                               Sfxiterator *csf,Env *env)
{
  env_error_check(env);
  while (true)
  {
    if (csf->phase == Partphase)
    {
      if (preparethispart(csf,mtime,env))
      {
        *len = csf->widthofpart;
        break;
      }
      if (hasspecialranges(csf->encseq))
      {
        csf->phase = Specialrangephase;
        csf->sri = newspecialrangeiterator(csf->encseq,
                                           ISDIRREVERSE(csf->readmode)
                                             ? false : true,
                                           env);
        csf->fusp.spaceSeqpos = csf->suftab;
        csf->fusp.allocatedSeqpos
          = CALLCASTFUNC(Seqpos,unsigned_long,
                         stpgetlargestwidth(csf->suftabparts));
        csf->fusp.nextfreeSeqpos = 0;
        csf->overhang.leftpos = csf->overhang.rightpos = 0;
      } else
      {
        csf->phase = Finalsuffixphase;
      }
    }
    if (csf->phase == Specialrangephase)
    {
      if (exhaustedspecialrangeiterator(csf->sri) && 
          csf->overhang.leftpos == csf->overhang.rightpos)
      {
        csf->phase = Finalsuffixphase;
      } else
      {
        fillspecialnextpage(csf);
        *len = (Seqpos) csf->fusp.nextfreeSeqpos;
        break;
      }
    }
  }
  return csf->suftab;
}
