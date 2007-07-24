/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <stdbool.h>
#ifndef NDEBUG
#define NDEBUG
#endif
#include <assert.h>
#include "libgtcore/env.h"
#include "libgtcore/strarray.h"
#include "chardef.h"
#include "spacedef.h"
#include "intcode-def.h"
#include "encseq-def.h"
#include "fbs-def.h"
#include "stamp.h"
#ifndef NDEBUG
#include "addnextchar.h"
#endif

#include "fbsadv.pr"

#include "readnextUchar.gen"

#ifdef SPECIALCASE4
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        if ((NUMOFCHARS) == DNAALPHASIZE)\
        {\
          CODE = MULT4((CODE) - MULTIMAPPOWER[(uint32_t) (LCHAR)]);\
        } else\
        {\
          CODE = ((CODE) - MULTIMAPPOWER[(uint32_t) (LCHAR)]) * (NUMOFCHARS);\
        }

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        if ((NUMOFCHARS) == DNAALPHASIZE)\
        {\
          CODE = MULT4((CODE) - MULTIMAPPOWER[(uint32_t) (LCHAR)]) | (CC);\
        } else\
        {\
          CODE = (Codetype) ((CODE) - MULTIMAPPOWER[(uint32_t) (LCHAR)]) *\
                            (NUMOFCHARS) + (CC);\
        }
#else
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        CODE = ((CODE) - MULTIMAPPOWER[(uint32_t) (LCHAR)]) * (NUMOFCHARS)

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        CODE = (Codetype) (((CODE) - MULTIMAPPOWER[(uint32_t) (LCHAR)]) * (NUMOFCHARS) + (CC))
#endif

#define ARRAY2DIMMALLOC(ARRAY2DIM, ROWS, COLUMNS, TYPE)\
        {\
          uint32_t rownumber;\
          ALLOCASSIGNSPACE(ARRAY2DIM,NULL,TYPE *,ROWS);\
          ALLOCASSIGNSPACE((ARRAY2DIM)[0],NULL,TYPE,(ROWS) * (COLUMNS));\
          for (rownumber = (uint32_t) 1; rownumber < (ROWS); rownumber++)\
          {\
            (ARRAY2DIM)[rownumber] = (ARRAY2DIM)[rownumber-1] + (COLUMNS);\
          }\
        }

#define ARRAY2DIMFREE(ARRAY2DIM)\
        FREESPACE((ARRAY2DIM)[0]);\
        FREESPACE(ARRAY2DIM)

#ifndef NDEBUG
static Codetype windowkmer2code(uint32_t numofchars,
                                uint32_t kmersize,
                                const Uchar *cyclicwindow,
                                uint32_t firstindex)
{
  uint32_t i;
  Codetype integercode;
  Uchar cc;
  bool foundspecial;

  cc = cyclicwindow[firstindex];
  if (ISSPECIAL(cc))
  {
    integercode = (Codetype) (numofchars-1);
    foundspecial = true;
  } else
  {
    integercode = (Codetype) cc;
    foundspecial = false;
  }
  for (i=(uint32_t) 1; i < kmersize; i++)
  {
    if (foundspecial)
    {
      ADDNEXTCHAR(integercode,numofchars-1,numofchars);
    } else
    {
      cc = cyclicwindow[(firstindex+i) % kmersize];
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
  return integercode;
}

static Codetype prefixwindowkmer2code(uint32_t firstspecialpos,
                                      uint32_t kmersize,
                                      Codetype **multimappower,
                                      const Uchar *cyclicwindow,
                                      uint32_t firstindex)
{
  uint32_t i;
  Codetype integercode = 0;
  Uchar cc;

  for (i=0; i<firstspecialpos; i++)
  {
    cc = cyclicwindow[(firstindex+i) % kmersize];
    integercode += multimappower[i][cc];
  }
  return integercode;
}

static Firstspecialpos determinefirstspecialposition(uint32_t windowwidth,
                                                     uint32_t kmersize,
                                                     const Uchar *cyclicwindow,
                                                     uint32_t firstindex)
{
  uint32_t i;
  Firstspecialpos fsp;

  for (i=0; i < windowwidth; i++)
  {
    if (ISSPECIAL(cyclicwindow[(firstindex+i) % kmersize]))
    {
      fsp.defined = true;
      fsp.specialpos = i;
      return fsp;
    }
  }
  fsp.defined = false;
  fsp.specialpos = 0; /* Just for satisfying the compiler */
  return fsp;
}
#endif

typedef struct
{
  uint32_t distvalue;
  Codetype codeforleftcontext;
} Queueelem;

typedef struct
{
  Queueelem *queuespace;  /* the space to store the queue elements */
  uint32_t enqueueindex,  /* entry into which element is to be enqued */
               dequeueindex,  /* last element of queue */
               queuesize,     /* size of the queue */
               noofelements;  /* no ofelements between enqueueindex+1 and 
                                 dequeindex */
} Specialpositions;

typedef struct
{
  Specialpositions spos;
  Uchar *cyclicwindow;
  uint32_t numofchars,
               kmersize,
               windowwidth,
               firstindex,
               *filltable,
               lengthwithoutspecial;
  Codetype codewithoutspecial,
           **multimappower;
} Streamstate;

static void specialemptyqueue(Specialpositions *spos,uint32_t queuesize,
                              Env *env)
{
  env_error_check(env);
  ALLOCASSIGNSPACE(spos->queuespace,NULL,Queueelem,queuesize);
  spos->noofelements = 0;
  spos->queuesize = queuesize;
  spos->dequeueindex = spos->enqueueindex = queuesize - 1;
}

static bool specialqueueisempty(const Specialpositions *spos)
{
  return (spos->noofelements == 0) ? true : false;
}

static Queueelem *specialheadofqueue(const Specialpositions *spos)
{
  return spos->queuespace + spos->dequeueindex;
}

static void specialdeleteheadofqueue(Specialpositions *spos)
{
  spos->noofelements--;
  if (spos->dequeueindex > 0)
  {
    spos->dequeueindex--;
  } else
  {
    spos->dequeueindex = spos->queuesize - 1;
  }
}

static void specialenqueue(Specialpositions *spos,Queueelem elem)
{
  spos->noofelements++;
  spos->queuespace[spos->enqueueindex] = elem;
  if (spos->enqueueindex > 0)
  {
    spos->enqueueindex--;
  } else
  {
    spos->enqueueindex = spos->queuesize - 1;
  }
}

static void specialwrapqueue(Specialpositions *spos,Env *env)
{
  env_error_check(env);
  FREESPACE(spos->queuespace);
}

static void updatespecialpositions(Streamstate *spwp,
                                   Uchar charcode,
                                   bool doshift,
                                   Uchar lchar)
{
  if (doshift)
  {
    if (!specialqueueisempty(&spwp->spos))
    {
      Queueelem *head;

      /* only here we add some element to the queue */
      head = specialheadofqueue(&spwp->spos);
      if (head->distvalue == 0)
      {
        specialdeleteheadofqueue(&spwp->spos);
        if (!specialqueueisempty(&spwp->spos))
        {
          head = specialheadofqueue(&spwp->spos);
          head->distvalue--;
        }
      } else
      {
        SUBTRACTLCHARANDSHIFT(head->codeforleftcontext,lchar,spwp->numofchars,
                              spwp->multimappower[0]);
        head->distvalue--;
      }
    }
  }
  if (ISSPECIAL(charcode))
  {
    /* only here we add some element to the queue */
    Queueelem newelem;

    if (specialqueueisempty(&spwp->spos))
    {
      newelem.distvalue = spwp->windowwidth-1;
    } else
    {
      newelem.distvalue = spwp->lengthwithoutspecial+1;
    }
    if (spwp->lengthwithoutspecial == spwp->kmersize)
    {
      SUBTRACTLCHARANDSHIFT(spwp->codewithoutspecial,lchar,
                            spwp->numofchars,spwp->multimappower[0]);
    }
    newelem.codeforleftcontext = spwp->codewithoutspecial;
    specialenqueue(&spwp->spos,newelem);
    spwp->lengthwithoutspecial = 0;
    spwp->codewithoutspecial = 0;
  } else
  {
    if (spwp->lengthwithoutspecial == spwp->kmersize)
    {
      SUBTRACTLCHARSHIFTADDNEXT(spwp->codewithoutspecial,
                                lchar,
                                spwp->numofchars,
                                spwp->multimappower[0],
                                charcode);
    } else
    {
      spwp->codewithoutspecial +=
        spwp->multimappower[spwp->lengthwithoutspecial][charcode];
      spwp->lengthwithoutspecial++;
    }
  }
}

static void shiftrightwithchar(
               void(*processkmercode)(void *,Codetype,Seqpos,
                                      const Firstspecialpos *,Env *),
               void *processkmercodeinfo,
               Streamstate *spwp,
               Seqpos currentposition,
               Uchar charcode,
               Env *env)
{
#ifndef NDEBUG
  Firstspecialpos firstspecialposbrute;
#endif

  env_error_check(env);
  if (spwp->windowwidth < spwp->kmersize)
  {
    spwp->windowwidth++;
    updatespecialpositions(spwp,charcode,false,0);
    spwp->cyclicwindow[spwp->windowwidth-1] = charcode;
  } else
  {
    updatespecialpositions(spwp,charcode,true,
                           spwp->cyclicwindow[spwp->firstindex]);
    spwp->cyclicwindow[spwp->firstindex] = charcode;
    if (spwp->firstindex == spwp->kmersize-1)
    {
      spwp->firstindex = 0;
    } else
    {
      spwp->firstindex++;
    }
  }
#ifndef NDEBUG
  if (!specialqueueisempty(&spwp->spos))
  {
    Queueelem *head = specialheadofqueue(&spwp->spos);
    Codetype tmpprefixcode = prefixwindowkmer2code(head->distvalue,
                                                   spwp->kmersize,
                                                   spwp->multimappower,
                                                   spwp->cyclicwindow,
                                                   spwp->firstindex);
    assert(tmpprefixcode == head->codeforleftcontext);
  }
  firstspecialposbrute = determinefirstspecialposition(spwp->windowwidth,
                                                       spwp->kmersize,
                                                       spwp->cyclicwindow,
                                                       spwp->firstindex);
  if (specialqueueisempty(&spwp->spos))
  {
    assert(!firstspecialposbrute.defined);
  } else
  {
    Queueelem *head = specialheadofqueue(&spwp->spos);
    assert(firstspecialposbrute.defined ? 1 : 0);
    assert(head->distvalue == firstspecialposbrute.specialpos);
  }
#endif
  if (spwp->windowwidth == spwp->kmersize)
  {
    Firstspecialpos localfirstspecial;
    Codetype code;

#ifndef NDEBUG
    Codetype wcode;

    wcode = windowkmer2code(spwp->numofchars,
                            spwp->kmersize,
                            spwp->cyclicwindow,
                            spwp->firstindex);
#endif
    if (specialqueueisempty(&spwp->spos))
    {
      localfirstspecial.defined = false;
      localfirstspecial.specialpos = 0;
      code = spwp->codewithoutspecial;
    } else
    {
      Queueelem *head = specialheadofqueue(&spwp->spos);
      code = head->codeforleftcontext + spwp->filltable[head->distvalue];
      localfirstspecial.defined = true;
      localfirstspecial.specialpos = head->distvalue;
    }
    assert(wcode == code);
    processkmercode(processkmercodeinfo,
                    code,
                    currentposition + 1 - spwp->kmersize,
                    &localfirstspecial,
                    env);
  }
}

static void initmultimappower(uint32_t ***multimappower,
                              uint32_t numofchars,
                              uint32_t kmersize,
                              Env *env)
{
  int offset;
  uint32_t thepower, mapindex, *mmptr;

  env_error_check(env);
  ARRAY2DIMMALLOC(*multimappower,kmersize,numofchars,uint32_t);
  thepower = (uint32_t) 1;
  for (offset=(int) (kmersize - 1); offset>=0; offset--)
  {
    mmptr = (*multimappower)[offset];
    mmptr[0] = 0;
    for (mapindex = (uint32_t) 1; mapindex < numofchars; mapindex++)
    {
      mmptr[mapindex] = mmptr[mapindex-1] + thepower;
    }
    thepower *= numofchars;
  }
}

static void filllargestchartable(uint32_t **filltable,
                                 uint32_t numofchars,
                                 uint32_t kmersize,
                                 Env *env)
{
  uint32_t *ptr;
  Codetype code;

  env_error_check(env);
  ALLOCASSIGNSPACE(*filltable,NULL,uint32_t,kmersize);
  code = numofchars;
  for (ptr = *filltable + kmersize - 1; ptr >= *filltable; ptr--)
  {
    *ptr = code-1;
    code *= numofchars;
  }
}

int getencseqkmersgeneric(
        const Encodedsequence *encseq,
        const StrArray *filenametab,
        void(*processkmercode)(void *,Codetype,Seqpos,
                               const Firstspecialpos *,Env *),
        void *processkmercodeinfo,
        uint32_t numofchars,
        uint32_t kmersize,
        const Uchar *symbolmap,
        bool plainformat,
        Env *env)
{
  uint32_t overshoot;
  Seqpos currentposition;
  Streamstate spwp;
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
  if(encseq != NULL)
  {
    Encodedsequencescanstate *esr;
    Seqpos totallength = getencseqtotallength(encseq);

    esr = initEncodedsequencescanstate(encseq,env);
    for (currentposition = 0; currentposition<totallength; currentposition++)
    {
      charcode = sequentialgetencodedchar(encseq,esr,currentposition);
      shiftrightwithchar(processkmercode,processkmercodeinfo,
                         &spwp,currentposition,charcode,env);
    }
    freeEncodedsequencescanstate(&esr,env);
  } else
  {
    Fastabufferstate fbs;
    int retval;

    initformatbufferstate(&fbs,
                          filenametab,
                          symbolmap,
                          plainformat,
                          NULL,
                          NULL,
                          env);
    for (currentposition = 0; /* Nothing */; currentposition++)
    {
      retval = readnextUchar(&charcode,&fbs,env);
      if (retval < 0)
      {
        return -1;
      }
      if (retval == 0)
      {
        break;
      }
      shiftrightwithchar(processkmercode,processkmercodeinfo,
                         &spwp,currentposition,charcode,env);
    }
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
  return 0;
}

int getfastastreamkmers(
        const StrArray *filenametab,
        void(*processkmercode)(void *,Codetype,Seqpos,
                               const Firstspecialpos *,Env *),
        void *processkmercodeinfo,
        uint32_t numofchars,
        uint32_t kmersize,
        const Uchar *symbolmap,
        bool plainformat,
        Env *env)
{
  return getencseqkmersgeneric(NULL,
                               filenametab,
                               processkmercode,
                               processkmercodeinfo,
                               numofchars,
                               kmersize,
                               symbolmap,
                               plainformat,
                               env);
}

void getencseqkmers(
        const Encodedsequence *encseq,
        void(*processkmercode)(void *,Codetype,Seqpos,
                               const Firstspecialpos *,Env *),
        void *processkmercodeinfo,
        uint32_t numofchars,
        uint32_t kmersize,
        Env *env)
{
  (void) getencseqkmersgeneric(encseq,
                               NULL,
                               processkmercode,
                               processkmercodeinfo,
                               numofchars,
                               kmersize,
                               NULL,
                               false,
                               env);
}
