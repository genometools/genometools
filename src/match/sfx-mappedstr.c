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
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/ma_api.h"
#include "core/checkencchar.h"
#include "core/encodedsequence_api.h"
#include "core/error_api.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/str_array.h"
#include "intcode-def.h"
#include "initbasepower.h"
#ifdef SKDEBUG
#include "sfx-nextchar.h"
#endif

#ifdef SPECIALCASE4
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        if ((NUMOFCHARS) == GT_DNAALPHASIZE)\
        {\
          CODE = MULT4((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]);\
        } else\
        {\
          CODE = ((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)])\
                 * (NUMOFCHARS);\
        }

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        if ((NUMOFCHARS) == GT_DNAALPHASIZE)\
        {\
          CODE = MULT4((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]) | (CC);\
        } else\
        {\
          CODE = (GtCodetype) ((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)])*\
                            (NUMOFCHARS) + (CC);\
        }
#else
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        CODE = ((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]) * (NUMOFCHARS)

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        CODE = (GtCodetype) (((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)])*\
                           (NUMOFCHARS) + (CC))
#endif

#ifdef SKDEBUG
static GtCodetype windowkmer2code(unsigned int numofchars,
                                unsigned int kmersize,
                                const GtUchar *cyclicwindow,
                                unsigned int firstindex)
{
  unsigned int i;
  GtCodetype integercode;
  GtUchar cc;
  bool foundspecial;

  cc = cyclicwindow[firstindex];
  if (ISSPECIAL(cc))
  {
    integercode = (GtCodetype) (numofchars-1);
    foundspecial = true;
  } else
  {
    integercode = (GtCodetype) cc;
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

static GtCodetype prefixwindowkmer2code(unsigned int firstspecialpos,
                                        unsigned int kmersize,
                                        const GtCodetype **multimappower,
                                        const GtUchar *cyclicwindow,
                                        unsigned int firstindex)
{
  unsigned int i;
  GtCodetype integercode = 0;
  GtUchar cc;

  for (i=0; i<firstspecialpos; i++)
  {
    cc = cyclicwindow[(firstindex+i) % kmersize];
    integercode += multimappower[i][cc];
  }
  return integercode;
}

static Firstspecialpos determinefirstspecialposition(unsigned int windowwidth,
                                                     unsigned int kmersize,
                                                     const GtUchar
                                                     *cyclicwindow,
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
  GtCodetype codeforleftcontext;
} Specialitem;

typedef struct
{
  Specialitem *queuespace;  /* the space to store the queue elements */
  unsigned int enqueueindex,  /* entry into which element is to be enqued */
               dequeueindex,  /* last element of queue */
               queuesize,     /* size of the queue */
               noofelements;  /* no ofelements between enqueueindex+1 and
                                 dequeindex */
} Specialpositions;

typedef struct
{
  Specialpositions spos;
  GtUchar cyclicwindow[MAXPREFIXLENGTH];
  unsigned int numofchars,
               kmersize,
               windowwidth,
               firstindex,
               lengthwithoutspecial;
  GtCodetype codewithoutspecial,
             *filltable,
             **multimappower;
} Kmerstream;

static void specialemptyqueue(Specialpositions *spos,unsigned int queuesize)
{
  spos->queuespace = gt_malloc(sizeof(*spos->queuespace) * queuesize);
  spos->noofelements = 0;
  spos->queuesize = queuesize;
  spos->dequeueindex = spos->enqueueindex = queuesize - 1;
}

static bool specialqueueisempty(const Specialpositions *spos)
{
  return (spos->noofelements == 0) ? true : false;
}

static Specialitem *specialheadofqueue(const Specialpositions *spos)
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

static void specialenqueue(Specialpositions *spos,Specialitem elem)
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
  gt_free(spos->queuespace);
}

static Kmerstream *kmerstream_new(unsigned int numofchars,unsigned int kmersize)
{
  Kmerstream *spwp;

  gt_assert(kmersize <= MAXPREFIXLENGTH);
  spwp = gt_malloc(sizeof (*spwp));
  spwp->multimappower = gt_initmultimappower(numofchars,kmersize);
  spwp->lengthwithoutspecial = 0;
  spwp->codewithoutspecial = 0;
  spwp->kmersize = kmersize;
  spwp->numofchars = numofchars;
  spwp->windowwidth = 0;
  spwp->firstindex = 0;
  specialemptyqueue(&spwp->spos,kmersize);
  spwp->filltable = gt_filllargestchartable(numofchars,kmersize);
  return spwp;
}

static void updatespecialpositions(Kmerstream *spwp,
                                   GtUchar charcode,
                                   bool doshift,
                                   GtUchar lchar)
{
  if (doshift)
  {
    if (!specialqueueisempty(&spwp->spos))
    {
      Specialitem *head;

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
    Specialitem newelem;

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
               void(*processkmercode)(void *,GtCodetype,unsigned long,
                                      const Firstspecialpos *),
               void *processkmercodeinfo,
               Kmerstream *spwp,
               unsigned long currentposition,
               GtUchar charcode)
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
    Specialitem *head = specialheadofqueue(&spwp->spos);
    GtCodetype tmpprefixcode = prefixwindowkmer2code(head->distvalue,
                                                   spwp->kmersize,
                                                   spwp->multimappower,
                                                   spwp->cyclicwindow,
                                                   spwp->firstindex);
    gt_assert(tmpprefixcode == head->codeforleftcontext);
  }
  firstspecialposbrute = determinefirstspecialposition(spwp->windowwidth,
                                                       spwp->kmersize,
                                                       spwp->cyclicwindow,
                                                       spwp->firstindex);
  if (specialqueueisempty(&spwp->spos))
  {
    gt_assert(!firstspecialposbrute.defined);
  } else
  {
    Specialitem *head = specialheadofqueue(&spwp->spos);
    gt_assert(firstspecialposbrute.defined ? 1 : 0);
    gt_assert(head->distvalue == firstspecialposbrute.specialpos);
  }
#endif
  if (spwp->windowwidth == spwp->kmersize)
  {
    Firstspecialpos localfirstspecial;
    GtCodetype code;

#ifdef SKDEBUG
    GtCodetype wcode;

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
      Specialitem *head = specialheadofqueue(&spwp->spos);
      code = head->codeforleftcontext + spwp->filltable[head->distvalue];
      localfirstspecial.defined = true;
      localfirstspecial.specialpos = head->distvalue;
    }
#ifdef SKDEBUG
    gt_assert(wcode == code);
#endif
    processkmercode(processkmercodeinfo,
                    code,
                    currentposition + 1 - spwp->kmersize,
                    &localfirstspecial);
  }
}

static void doovershoot(Kmerstream *spwp,
                        void(*processkmercode)(void *,GtCodetype,unsigned long,
                                               const Firstspecialpos *),
                        void *processkmercodeinfo,
                        unsigned long currentposition,
                        unsigned int kmersize)
{
  unsigned int overshoot;

  for (overshoot=0; overshoot<kmersize; overshoot++)
  {
    shiftrightwithchar(processkmercode,processkmercodeinfo,spwp,
                       currentposition + overshoot,(GtUchar) WILDCARD);
  }
}

static void freestreamstate(Kmerstream *spwp)
{
  gt_free(spwp->filltable);
  gt_multimappowerfree(&spwp->multimappower);
  specialwrapqueue(&spwp->spos);
  gt_free(spwp);
}

void getencseqkmers(
        const GtEncodedsequence *encseq,
        GtReadmode readmode,
        void(*processkmercode)(void *,
                               GtCodetype,
                               unsigned long,
                               const Firstspecialpos *),
        void *processkmercodeinfo,
        unsigned int kmersize)
{
  unsigned long currentposition = 0, totallength;
  Kmerstream *spwp;
  GtUchar charcode;
  GtEncodedsequenceScanstate *esr;
  unsigned int numofchars;

  totallength = gt_encodedsequence_total_length(encseq);
  numofchars = gt_alphabet_num_of_chars(gt_encodedsequence_alphabet(encseq));
  spwp = kmerstream_new(numofchars,kmersize);
  esr = gt_encodedsequence_scanstate_new(encseq,readmode,0);
  for (currentposition = 0; currentposition<totallength; currentposition++)
  {
    charcode = gt_encodedsequence_get_encoded_char_sequential(encseq,esr,
                                                              currentposition,
                                                              readmode);
    GT_CHECKENCCHAR(charcode,encseq,currentposition,readmode);
    shiftrightwithchar(processkmercode,processkmercodeinfo,
                       spwp,currentposition,charcode);
  }
  gt_encodedsequence_scanstate_delete(esr);
  doovershoot(spwp,
              processkmercode,
              processkmercodeinfo,
              currentposition,
              kmersize);
  freestreamstate(spwp);
}

int getfastastreamkmers(
        const GtStrArray *filenametab,
        void(*processkmercode)(void *,
                               GtCodetype,
                               unsigned long,
                               const Firstspecialpos *),
        void *processkmercodeinfo,
        unsigned int numofchars,
        unsigned int kmersize,
        const GtUchar *symbolmap,
        bool plainformat,
        GtError *err)
{
  unsigned long currentposition = 0;
  Kmerstream *spwp;
  GtUchar charcode;
  bool haserr = false;
  GtSequenceBuffer *fb;
  int retval;

  gt_error_check(err);
  spwp = kmerstream_new(numofchars,kmersize);
  if (plainformat)
  {
    fb = gt_sequence_buffer_plain_new(filenametab);
  } else
  {
    fb = gt_sequence_buffer_new_guess_type(filenametab, err);
  }
  if (fb == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    gt_sequence_buffer_set_symbolmap(fb, symbolmap);

    for (currentposition = 0; /* Nothing */; currentposition++)
    {
      retval = gt_sequence_buffer_next(fb,&charcode,err);
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
                         spwp,currentposition,charcode);
    }
    gt_sequence_buffer_delete(fb);
  }
  if (!haserr)
  {
    doovershoot(spwp,
                processkmercode,
                processkmercodeinfo,
                currentposition,
                kmersize);
  }
  freestreamstate(spwp);
  return haserr ? -1 : 0;
}
