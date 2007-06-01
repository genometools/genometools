#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <stdbool.h>
#define NDEBUG
#include <assert.h>
#include "libgtcore/env.h"
#include "types.h"
#include "chardef.h"
#include "inputsymbol.h"
#include "spacedef.h"
#include "genstream.h"
#include "encseq-def.h"
#ifndef NDEBUG
#include "addnextchar.h"
#endif

#include "genericstream.pr"

#ifdef SPECIALCASE4
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        if ((NUMOFCHARS) == DNAALPHASIZE)\
        {\
          CODE = MULT4((CODE) - MULTIMAPPOWER[(Uint) (LCHAR)]);\
        } else\
        {\
          CODE = ((CODE) - MULTIMAPPOWER[(Uint) (LCHAR)]) * (NUMOFCHARS);\
        }

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        if ((NUMOFCHARS) == DNAALPHASIZE)\
        {\
          CODE = MULT4((CODE) - MULTIMAPPOWER[(Uint) (LCHAR)]) | (CC);\
        } else\
        {\
          CODE = ((CODE) - MULTIMAPPOWER[(Uint) (LCHAR)]) * (NUMOFCHARS)\
                  + (CC);\
        }
#else
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        CODE = ((CODE) - MULTIMAPPOWER[(Uint) (LCHAR)]) * (NUMOFCHARS)

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        CODE = ((CODE) - MULTIMAPPOWER[(Uint) (LCHAR)]) * (NUMOFCHARS) + (CC)
#endif

#define ARRAY2DIMMALLOC(ARRAY2DIM, ROWS, COLUMNS, TYPE)\
        {\
          Uint rownumber;\
          ALLOCASSIGNSPACE(ARRAY2DIM,NULL,TYPE *,ROWS);\
          ALLOCASSIGNSPACE((ARRAY2DIM)[0],NULL,TYPE,(ROWS) * (COLUMNS));\
          for (rownumber = UintConst(1); rownumber < (ROWS); rownumber++)\
          {\
            (ARRAY2DIM)[rownumber] = (ARRAY2DIM)[rownumber-1] + (COLUMNS);\
          }\
        }

#define ARRAY2DIMFREE(ARRAY2DIM)\
        FREESPACE((ARRAY2DIM)[0]);\
        FREESPACE(ARRAY2DIM)

#ifndef NDEBUG
static Uint windowkmer2code(Uint numofchars,
                            Uint kmersize,
                            const Uchar *cyclicwindow,
                            Uint firstindex)
{
  Uint i, integercode;
  Uchar cc;
  bool foundspecial;

  cc = cyclicwindow[firstindex];
  if (ISSPECIAL(cc))
  {
    integercode = (Uint) (numofchars-1);
    foundspecial = true;
  } else
  {
    integercode = (Uint) cc;
    foundspecial = false;
  }
  for (i=UintConst(1); i < kmersize; i++)
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

static Uint prefixwindowkmer2code(Uint firstspecialpos,
                                  Uint kmersize,
                                  Uint **multimappower,
                                  const Uchar *cyclicwindow,
                                  Uint firstindex)
{
  Uint i, integercode = 0;
  Uchar cc;

  for (i=0; i<firstspecialpos; i++)
  {
    cc = cyclicwindow[(firstindex+i) % kmersize];
    integercode += multimappower[i][cc];
  }
  return integercode;
}

static DefinedUint determinefirstspecialposition(Uint windowwidth,
                                                 Uint kmersize,
                                                 const Uchar *cyclicwindow,
                                                 Uint firstindex)
{
  Uint i;
  DefinedUint retval;

  for (i=0; i < windowwidth; i++)
  {
    if (ISSPECIAL(cyclicwindow[(firstindex+i) % kmersize]))
    {
      retval.defined = true;
      retval.uintvalue = i;
      return retval;
    }
  }
  retval.defined = false;
  retval.uintvalue = 0;
  return retval;
}
#endif

typedef struct
{
  Uint distvalue,
       codeforleftcontext;
} Queueelem;

typedef struct
{
  Queueelem *queuespace;  /* the space to store the queue elements */
  Uint enqueueindex,  /* points to entry into which element is to be enqued */
       dequeueindex,  /* last element of queue */
       queuesize,     /* size of the queue */
       noofelements;  /* no ofelements between enqueueindex+1 and dequeindex */
} Specialpositions;

typedef struct
{
  Specialpositions spos;
  Uchar *cyclicwindow;
  Uint numofchars,
       kmersize,
       windowwidth,
       firstindex,
       lengthwithoutspecial,
       codewithoutspecial,
       *filltable,
       **multimappower;
} Streamstate;

static void specialemptyqueue(Specialpositions *spos,Uint queuesize, Env *env)
{
  env_error_check(env);
  ALLOCASSIGNSPACE(spos->queuespace,NULL,Queueelem,queuesize);
  spos->noofelements = 0;
  spos->queuesize = queuesize;
  spos->dequeueindex
    = spos->enqueueindex = queuesize - 1;
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
      SUBTRACTLCHARSHIFTADDNEXT(spwp->codewithoutspecial,lchar,
                                spwp->numofchars,spwp->multimappower[0],
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
               void(*processkmercode)(void *,Uint,Uint64,
                                      const DefinedUint *,Env *),
               void *processkmercodeinfo,
               Streamstate *spwp,
               Uint64 currentposition,
               Uchar charcode,
               Env *env)
{
#ifndef NDEBUG
  DefinedUint firstspecialposbrute;
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
    Uint tmpprefixcode = prefixwindowkmer2code(head->distvalue,
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
    assert(head->distvalue == firstspecialposbrute.uintvalue);
  }
#endif
  if (spwp->windowwidth == spwp->kmersize)
  {
    DefinedUint localfirstspecial;
    Uint code;

#ifndef NDEBUG
    Uint wcode;

    wcode = windowkmer2code(spwp->numofchars,
                            spwp->kmersize,
                            spwp->cyclicwindow,
                            spwp->firstindex);
#endif
    if (specialqueueisempty(&spwp->spos))
    {
      localfirstspecial.defined = false;
      localfirstspecial.uintvalue = 0;
      code = spwp->codewithoutspecial;
    } else
    {
      Queueelem *head = specialheadofqueue(&spwp->spos);
      code = head->codeforleftcontext + spwp->filltable[head->distvalue];
      localfirstspecial.defined = true;
      localfirstspecial.uintvalue = head->distvalue;
    }
    assert(wcode == code);
    processkmercode(processkmercodeinfo,
                    code,
                    currentposition + 1 - spwp->kmersize,
                    &localfirstspecial,
                    env);
  }
}

static void initmultimappower(Uint ***multimappower,Uint numofchars,
                              unsigned int kmersize,Env *env)
{
  int offset;
  Uint thepower, mapindex, *mmptr;

  env_error_check(env);
  ARRAY2DIMMALLOC(*multimappower,kmersize,numofchars,Uint);
  thepower = UintConst(1);
  for (offset=(int) (kmersize - 1); offset>=0; offset--)
  {
    mmptr = (*multimappower)[offset];
    mmptr[0] = 0;
    for (mapindex = UintConst(1); mapindex < numofchars; mapindex++)
    {
      mmptr[mapindex] = mmptr[mapindex-1] + thepower;
    }
    thepower *= numofchars;
  }
}

static void filllargestchartable(Uint **filltable,Uint numofchars,
                                 unsigned int kmersize,Env *env)
{
  Uint *ptr, code;

  env_error_check(env);
  ALLOCASSIGNSPACE(*filltable,NULL,Uint,kmersize);
  code = numofchars;
  for (ptr = *filltable + kmersize - 1; ptr >= *filltable; ptr--)
  {
    *ptr = code-1;
    code *= numofchars;
  }
}

int getfastastreamkmers(
        const char **filenametab,
        unsigned int numoffiles,
        void(*processkmercode)(void *,Uint,Uint64,const DefinedUint *,Env *),
        void *processkmercodeinfo,
        Uint numofchars,
        unsigned int kmersize,
        const Uchar *symbolmap,
        Env *env)
{
  unsigned int filenum;
  Fgetcreturntype currentchar;
  bool indesc, firstseq = true;
  unsigned int overshoot;
  Uint linenum = UintConst(1);
  Uint64 currentposition = 0;
  Streamstate spwp;
  Genericstream inputstream;
  Uchar charcode;

  env_error_check(env);
  assert(numoffiles > 0);
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
  for (filenum = 0; filenum < numoffiles; filenum++)
  {
    opengenericstream(&inputstream,filenametab[filenum]);
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
              /* first separator symbol */
              shiftrightwithchar(processkmercode,processkmercodeinfo,
                                 &spwp,currentposition,(Uchar) SEPARATOR,env);
              currentposition++;
            }
            indesc = true;
          } else
          {
            charcode = symbolmap[(Uint) currentchar];
            if (charcode == (Uchar) UNDEFCHAR)
            {
              env_error_set(env,
                            "illegal character '%c': file \"%s\", line %lu",
                            currentchar,
                            filenametab[filenum],
                            (Showuint) linenum);
              return -1;
            }
            shiftrightwithchar(processkmercode,processkmercodeinfo,
                               &spwp,currentposition,charcode,env);
            currentposition++;
          }
        }
      }
    }
    closegenericstream(&inputstream,filenametab[filenum]);
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
                  filenametab[0]);
    return -2;
  }
  return 0;
}

void getencseqkmers(
        const Encodedsequence *encseq,
        Uint64 totallength,
        void(*processkmercode)(void *,Uint,Uint64,const DefinedUint *,Env *),
        void *processkmercodeinfo,
        Uint numofchars,
        unsigned int kmersize,
        Env *env)
{
  unsigned int overshoot;
  Uint64 currentposition;
  Streamstate spwp;
  Uchar charcode;
  Encodedsequencescanstate *esr;

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
  esr = initEncodedsequencescanstate(encseq,env);
  for (currentposition = 0; currentposition<totallength; currentposition++)
  {
    charcode = sequentialgetencodedchar64(encseq,esr,currentposition);
    shiftrightwithchar(processkmercode,processkmercodeinfo,
                       &spwp,currentposition,charcode,env);
  }
  freeEncodedsequencescanstate(&esr,env);
  for (overshoot=0; overshoot<kmersize; overshoot++)
  {
    shiftrightwithchar(processkmercode,processkmercodeinfo,&spwp,
                       currentposition + overshoot,(Uchar) WILDCARD,env);
  }
  FREESPACE(spwp.cyclicwindow);
  FREESPACE(spwp.filltable);
  ARRAY2DIMFREE(spwp.multimappower);
  specialwrapqueue(&spwp.spos,env);
}
