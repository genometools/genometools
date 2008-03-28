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
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <stdbool.h>
#include <assert.h>
#include "libgtcore/chardef.h"
#include "libgtcore/error.h"
#include "libgtcore/fastabuffer.h"
#include "libgtcore/strarray.h"
#include "spacedef.h"
#include "intcode-def.h"
#include "encseq-def.h"
#include "stamp.h"
#ifdef SKDEBUG
#include "sfx-nextchar.h"
#endif

#include "initbasepower.pr"

#ifdef SPECIALCASE4
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        if ((NUMOFCHARS) == DNAALPHASIZE)\
        {\
          CODE = MULT4((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]);\
        } else\
        {\
          CODE = ((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)])\
                 * (NUMOFCHARS);\
        }

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        if ((NUMOFCHARS) == DNAALPHASIZE)\
        {\
          CODE = MULT4((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]) | (CC);\
        } else\
        {\
          CODE = (Codetype) ((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]) *\
                            (NUMOFCHARS) + (CC);\
        }
#else
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        CODE = ((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]) * (NUMOFCHARS)

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        CODE = (Codetype) (((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]) *\
                           (NUMOFCHARS) + (CC))
#endif

#ifdef SKDEBUG
static Codetype windowkmer2code(unsigned int numofchars,
                                unsigned int kmersize,
                                const Uchar *cyclicwindow,
                                unsigned int firstindex)
{
  unsigned int i;
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
  for (i=1U; i < kmersize; i++)
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

static Codetype prefixwindowkmer2code(unsigned int firstspecialpos,
                                      unsigned int kmersize,
                                      const Codetype **multimappower,
                                      const Uchar *cyclicwindow,
                                      unsigned int firstindex)
{
  unsigned int i;
  Codetype integercode = 0;
  Uchar cc;

  for (i=0; i<firstspecialpos; i++)
  {
    cc = cyclicwindow[(firstindex+i) % kmersize];
    integercode += multimappower[i][cc];
  }
  return integercode;
}

static Firstspecialpos determinefirstspecialposition(unsigned int windowwidth,
                                                     unsigned int kmersize,
                                                     const Uchar *cyclicwindow,
                                                     unsigned int firstindex)
{
  unsigned int i;
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
  unsigned int distvalue;
  Codetype codeforleftcontext;
} Queueelem;

typedef struct
{
  Queueelem *queuespace;  /* the space to store the queue elements */
  unsigned int enqueueindex,  /* entry into which element is to be enqued */
               dequeueindex,  /* last element of queue */
               queuesize,     /* size of the queue */
               noofelements;  /* no ofelements between enqueueindex+1 and
                                 dequeindex */
} Specialpositions;

typedef struct
{
  Specialpositions spos;
  Uchar *cyclicwindow;
  unsigned int numofchars,
               kmersize,
               windowwidth,
               firstindex,
               lengthwithoutspecial;
  Codetype codewithoutspecial,
           *filltable;
  Codetype **multimappower;
} Streamstate;

static void specialemptyqueue(Specialpositions *spos,unsigned int queuesize)
{
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

static void specialwrapqueue(Specialpositions *spos)
{
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
                                      const Firstspecialpos *),
               void *processkmercodeinfo,
               Streamstate *spwp,
               Seqpos currentposition,
               Uchar charcode)
{
#ifdef SKDEBUG
  Firstspecialpos firstspecialposbrute;
#endif

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
#ifdef SKDEBUG
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

#ifdef SKDEBUG
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
#ifdef SKDEBUG
    assert(wcode == code);
#endif
    processkmercode(processkmercodeinfo,
                    code,
                    currentposition + 1 - spwp->kmersize,
                    &localfirstspecial);
  }
}

static void filllargestchartable(Codetype **filltable,
                                 unsigned int numofchars,
                                 unsigned int kmersize)
{
  Codetype code, *ptr;

  ALLOCASSIGNSPACE(*filltable,NULL,Codetype,kmersize);
  code = (Codetype) numofchars;
  for (ptr = *filltable + kmersize - 1; ptr >= *filltable; ptr--)
  {
    *ptr = code-1;
    code *= numofchars;
  }
}

static void initstreamstate(Streamstate *spwp,unsigned int numofchars,
                            unsigned int kmersize)
{
  spwp->multimappower = initmultimappower(numofchars,kmersize);
  spwp->lengthwithoutspecial = 0;
  spwp->codewithoutspecial = 0;
  spwp->kmersize = kmersize;
  spwp->numofchars = numofchars;
  spwp->windowwidth = 0;
  spwp->firstindex = 0;
  ALLOCASSIGNSPACE(spwp->cyclicwindow,NULL,Uchar,kmersize);
  specialemptyqueue(&spwp->spos,kmersize);
  filllargestchartable(&spwp->filltable,numofchars,kmersize);
}

static void freestreamstate(Streamstate *spwp)
{
  FREESPACE(spwp->cyclicwindow);
  FREESPACE(spwp->filltable);
  multimappowerfree(&spwp->multimappower);
  specialwrapqueue(&spwp->spos);
}

static void doovershoot(Streamstate *spwp,
                        void(*processkmercode)(void *,Codetype,Seqpos,
                                               const Firstspecialpos *),
                        void *processkmercodeinfo,
                        Seqpos currentposition,
                        unsigned int kmersize)
{
  unsigned int overshoot;

  for (overshoot=0; overshoot<kmersize; overshoot++)
  {
    shiftrightwithchar(processkmercode,processkmercodeinfo,spwp,
                       currentposition + overshoot,(Uchar) WILDCARD);
  }
}

void getencseqkmers(
        const Encodedsequence *encseq,
        Readmode readmode,
        void(*processkmercode)(void *,Codetype,Seqpos,const Firstspecialpos *),
        void *processkmercodeinfo,
        unsigned int numofchars,
        unsigned int kmersize)
{
  Seqpos currentposition = 0, totallength;
  Streamstate spwp;
  Uchar charcode;
  Encodedsequencescanstate *esr;

  totallength = getencseqtotallength(encseq);
  initstreamstate(&spwp,numofchars,kmersize);
  esr = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr,encseq,readmode,0);
  for (currentposition = 0; currentposition<totallength; currentposition++)
  {
    charcode = sequentialgetencodedchar(encseq,esr,currentposition,readmode);
    CHECKENCCHAR(charcode,encseq,currentposition,readmode);
    shiftrightwithchar(processkmercode,processkmercodeinfo,
                       &spwp,currentposition,charcode);
  }
  if (esr != NULL)
  {
    freeEncodedsequencescanstate(&esr);
  }
  doovershoot(&spwp,
              processkmercode,
              processkmercodeinfo,
              currentposition,
              kmersize);
  freestreamstate(&spwp);
}

int getfastastreamkmers(
        const StrArray *filenametab,
        void(*processkmercode)(void *,Codetype,Seqpos,const Firstspecialpos *),
        void *processkmercodeinfo,
        unsigned int numofchars,
        unsigned int kmersize,
        const Uchar *symbolmap,
        bool plainformat,
        Error *err)
{
  Seqpos currentposition = 0;
  Streamstate spwp;
  Uchar charcode;
  bool haserr = false;
  FastaBuffer *fb;
  int retval;

  error_check(err);
  initstreamstate(&spwp,numofchars,kmersize);
  fb = fastabuffer_new(filenametab,
                       symbolmap,
                       plainformat,
                       NULL,
                       NULL,
                       NULL);
  for (currentposition = 0; /* Nothing */; currentposition++)
  {
    retval = fastabuffer_next(fb,&charcode,err);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    shiftrightwithchar(processkmercode,processkmercodeinfo,
                       &spwp,currentposition,charcode);
  }
  fastabuffer_delete(fb);
  if (!haserr)
  {
    doovershoot(&spwp,
                processkmercode,
                processkmercodeinfo,
                currentposition,
                kmersize);
  }
  freestreamstate(&spwp);
  return haserr ? -1 : 0;
}
