/*
  Copyright (c) 2007-2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef S_SPLINT_S
#include <ctype.h>
#include <zlib.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/ma_api.h"
#include "core/checkencchar.h"
#include "core/divmodmul.h"
#include "core/alphabet.h"
#include "core/error_api.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/str_array.h"
#include "core/intbits.h"
#include "core/encseq.h"
#include "core/format64.h"
#include "intcode-def.h"
#include "initbasepower.h"
#include "sfx-mappedstr.h"
#undef SKDEBUG
#ifdef SKDEBUG
#include "sfx-nextchar.h"
#endif

#undef SPECIALCASE4
#ifdef SPECIALCASE4
#define SUBTRACTLCHARANDSHIFT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER)\
        if ((NUMOFCHARS) == GT_DNAALPHASIZE)\
        {\
          CODE = GT_MULT4((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]);\
        } else\
        {\
          CODE = ((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)])\
                 * (NUMOFCHARS);\
        }

#define SUBTRACTLCHARSHIFTADDNEXT(CODE,LCHAR,NUMOFCHARS,MULTIMAPPOWER,CC)\
        if ((NUMOFCHARS) == GT_DNAALPHASIZE)\
        {\
          CODE = GT_MULT4((CODE) - MULTIMAPPOWER[(unsigned int) (LCHAR)]) |\
                 (CC);\
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

static bool determinefirstspecialposition(unsigned int *firstspecialpos,
                                          unsigned int windowwidth,
                                          unsigned int kmersize,
                                                       const GtUchar
                                                         *cyclicwindow,
                                                       unsigned int firstindex)
{
  unsigned int i;

  for (i=0; i < windowwidth; i++)
  {
    if (ISSPECIAL(cyclicwindow[(firstindex+i) % kmersize]))
    {
      *firstspecialpos = i;
      return true;
    }
  }
  *firstspecialpos = 0; /* Just for satisfying the compiler */
  return false;
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
  GtKmercode currentkmercode;
} Kmerstream;

static void specialemptyqueue(Specialpositions *spos,unsigned int queuesize)
{
  spos->queuespace = gt_malloc(sizeof (*spos->queuespace) * queuesize);
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

  gt_assert(kmersize <= (unsigned int) MAXPREFIXLENGTH);
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
      if (head->distvalue > 0)
      {
        SUBTRACTLCHARANDSHIFT(head->codeforleftcontext,lchar,spwp->numofchars,
                              spwp->multimappower[0]);
        head->distvalue--;
      } else
      {
        specialdeleteheadofqueue(&spwp->spos);
        if (!specialqueueisempty(&spwp->spos))
        {
          head = specialheadofqueue(&spwp->spos);
          head->distvalue--;
        }
      }
    }
  }
  if (ISNOTSPECIAL(charcode))
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
  } else
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
  }
}

static void kmerstream_newcode(GtKmercode *kmercode, Kmerstream *spwp)
{
#ifdef SKDEBUG
  bool firstspecialbrutedefined;
  unsigned int firstspecialbrute;

  if (!specialqueueisempty(&spwp->spos))
  {
    Specialitem *head = specialheadofqueue(&spwp->spos);
    GtCodetype tmpprefixcode = prefixwindowkmer2code(head->distvalue,
                                                     spwp->kmersize,
                                                     (const GtCodetype **)
                                                       spwp->multimappower,
                                                     spwp->cyclicwindow,
                                                     spwp->firstindex);
    gt_assert(tmpprefixcode == head->codeforleftcontext);
  }
  firstspecialbrutedefined
    = determinefirstspecialposition(&firstspecialbrute,
                                    spwp->windowwidth,
                                    spwp->kmersize,
                                    spwp->cyclicwindow,
                                    spwp->firstindex);
  if (specialqueueisempty(&spwp->spos))
  {
    gt_assert(!firstspecialbrutedefined);
  } else
  {
    Specialitem *head = specialheadofqueue(&spwp->spos);
    gt_assert(firstspecialbrutedefined ? 1 : 0);
    gt_assert(head->distvalue == firstspecialbrute);
  }
#endif
  {
#ifdef SKDEBUG
    GtCodetype wcode;

    wcode = windowkmer2code(spwp->numofchars,
                            spwp->kmersize,
                            spwp->cyclicwindow,
                            spwp->firstindex);
#endif
    if (specialqueueisempty(&spwp->spos))
    {
      kmercode->definedspecialposition = false;
      kmercode->specialposition = 0;
      kmercode->code = spwp->codewithoutspecial;
    } else
    {
      Specialitem *head = specialheadofqueue(&spwp->spos);
      kmercode->code = head->codeforleftcontext +
                       spwp->filltable[head->distvalue];
      kmercode->definedspecialposition = true;
      kmercode->specialposition = head->distvalue;
    }
#ifdef SKDEBUG
    gt_assert(wcode == kmercode->code);
#endif
  }
}

static void shiftrightwithchar(Kmerstream *spwp,GtUchar charcode)
{
  gt_assert(spwp->windowwidth == spwp->kmersize);
  updatespecialpositions(spwp,charcode,true,
                         spwp->cyclicwindow[spwp->firstindex]);
  spwp->cyclicwindow[spwp->firstindex] = charcode;
  if (spwp->firstindex < spwp->kmersize-1)
  {
    spwp->firstindex++;
  } else
  {
    spwp->firstindex = 0;
  }
}

static void kmerstream_delete(Kmerstream *spwp)
{
  gt_free(spwp->filltable);
  gt_multimappower_delete(spwp->multimappower);
  specialwrapqueue(&spwp->spos);
  gt_free(spwp);
}

struct GtKmercodeiterator
{
  unsigned long totallength, startpos;
  const GtEncseq *encseq;
  GtEncseqReader *esr;
  GtReadmode readmode;
  Kmerstream *spwp;
  unsigned long currentposition;
  bool hasprocessedfirst, inputexhausted;
  GtKmercode kmercode;
  GtSequenceBuffer *fb; /* only for generating from file */
};

/*@notnull@*/ GtKmercodeiterator *gt_kmercodeiterator_encseq_new(
                                            const GtEncseq *encseq,
                                            GtReadmode readmode,
                                            unsigned int kmersize,
                                            unsigned long startpos)
{
  GtKmercodeiterator *kmercodeiterator;
  unsigned int numofchars;
  GtUchar charcode;

  gt_assert(!GT_ISDIRREVERSE(readmode) || startpos == 0);
  kmercodeiterator = gt_malloc(sizeof (*kmercodeiterator));
  kmercodeiterator->totallength = gt_encseq_total_length(encseq);
  kmercodeiterator->startpos = startpos;
  gt_assert(startpos < kmercodeiterator->totallength);
  if (kmercodeiterator->totallength - startpos < (unsigned long) kmersize)
  {
    kmercodeiterator->inputexhausted = true;
    kmercodeiterator->fb = NULL;
    kmercodeiterator->encseq = encseq;
    kmercodeiterator->esr = NULL;
    kmercodeiterator->spwp = NULL;
  } else
  {
    kmercodeiterator->inputexhausted = false;
    kmercodeiterator->fb = NULL;
    kmercodeiterator->encseq = encseq;
    kmercodeiterator->readmode = readmode;
    kmercodeiterator->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                                  readmode,
                                                                  startpos);
    numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
    kmercodeiterator->spwp = kmerstream_new(numofchars,kmersize);
    kmercodeiterator->hasprocessedfirst = false;
    for (kmercodeiterator->currentposition = startpos;
         kmercodeiterator->currentposition < startpos+(unsigned long) kmersize;
         kmercodeiterator->currentposition++)
    {
      charcode = gt_encseq_reader_next_encoded_char(kmercodeiterator->esr);
      kmercodeiterator->spwp->windowwidth++;
      updatespecialpositions(kmercodeiterator->spwp,charcode,false,0);
      kmercodeiterator->spwp->cyclicwindow[kmercodeiterator->
                                           spwp->windowwidth-1] = charcode;
    }
  }
  return kmercodeiterator;
}

const GtKmercode *gt_kmercodeiterator_encseq_next(
                             GtKmercodeiterator *kmercodeiterator)
{
  if (!kmercodeiterator->hasprocessedfirst)
  {
    gt_assert(kmercodeiterator->currentposition
              == kmercodeiterator->startpos +
                   (unsigned long) kmercodeiterator->spwp->kmersize);
    kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
    kmercodeiterator->hasprocessedfirst = true;
    return &kmercodeiterator->kmercode;
  }
  if (kmercodeiterator->currentposition < kmercodeiterator->totallength)
  {
    GtUchar charcode
      = gt_encseq_reader_next_encoded_char(kmercodeiterator->esr);
    shiftrightwithchar(kmercodeiterator->spwp,charcode);
    kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
    kmercodeiterator->currentposition++;
    return &kmercodeiterator->kmercode;
  }
  if (kmercodeiterator->currentposition < kmercodeiterator->totallength +
                                          kmercodeiterator->spwp->kmersize)
  {
    shiftrightwithchar(kmercodeiterator->spwp,(GtUchar) WILDCARD);
    kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
    kmercodeiterator->currentposition++;
    return &kmercodeiterator->kmercode;
  }
  return NULL;
}

const GtKmercode *gt_kmercodeiterator_encseq_nonspecial_next(
                             GtKmercodeiterator *kmercodeiterator)
{
  while (true)
  {
    if (!kmercodeiterator->hasprocessedfirst)
    {
      gt_assert(kmercodeiterator->currentposition
                == kmercodeiterator->startpos +
                     (unsigned long) kmercodeiterator->spwp->kmersize);
      kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
      kmercodeiterator->hasprocessedfirst = true;
      if (!kmercodeiterator->kmercode.definedspecialposition)
      {
        return &kmercodeiterator->kmercode;
      }
    } else
    {
      if (kmercodeiterator->currentposition < kmercodeiterator->totallength)
      {
        GtUchar charcode
          = gt_encseq_reader_next_encoded_char(kmercodeiterator->esr);
        shiftrightwithchar(kmercodeiterator->spwp,charcode);
        kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
        kmercodeiterator->currentposition++;
        if (!kmercodeiterator->kmercode.definedspecialposition)
        {
          return &kmercodeiterator->kmercode;
        }
      } else
      {
        break;
      }
    }
  }
  return NULL;
}

GtKmercodeiterator *gt_kmercodeiterator_filetab_new(
                                                const GtStrArray *filenametab,
                                                unsigned int numofchars,
                                                unsigned int kmersize,
                                                const GtUchar *symbolmap,
                                                bool plainformat,
                                                GtError *err)
{
  GtKmercodeiterator *kmercodeiterator;
  GtUchar charcode;
  bool haserr = false;
  int retval;

  gt_error_check(err);
  kmercodeiterator = gt_malloc(sizeof (*kmercodeiterator));
  kmercodeiterator->esr = NULL;
  kmercodeiterator->hasprocessedfirst = false;
  kmercodeiterator->inputexhausted = false;
  kmercodeiterator->spwp = kmerstream_new(numofchars,kmersize);
  kmercodeiterator->totallength = 0;
  if (plainformat)
  {
    kmercodeiterator->fb = gt_sequence_buffer_plain_new(filenametab);
  } else
  {
    kmercodeiterator->fb = gt_sequence_buffer_new_guess_type(filenametab, err);
  }
  if (kmercodeiterator->fb == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    gt_sequence_buffer_set_symbolmap(kmercodeiterator->fb, symbolmap);
    for (kmercodeiterator->currentposition = 0;
         kmercodeiterator->currentposition < (unsigned long) kmersize;
         kmercodeiterator->currentposition++)
    {
      retval = gt_sequence_buffer_next(kmercodeiterator->fb,&charcode,err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        kmercodeiterator->inputexhausted = true;
        break;
      }
      kmercodeiterator->spwp->windowwidth++;
      updatespecialpositions(kmercodeiterator->spwp,charcode,false,0);
      kmercodeiterator->spwp->cyclicwindow[kmercodeiterator->
                                           spwp->windowwidth-1] = charcode;
    }
  }
  if (haserr)
  {
    gt_kmercodeiterator_delete(kmercodeiterator);
    return NULL;
  }
  return kmercodeiterator;
}

int gt_kmercodeiterator_filetab_next(const GtKmercode **kmercodeptr,
                                     GtKmercodeiterator *kmercodeiterator,
                                     GtError *err)
{
  if (!kmercodeiterator->inputexhausted)
  {
    if (kmercodeiterator->hasprocessedfirst)
    {
      GtUchar charcode;
      int retval;

      retval = gt_sequence_buffer_next(kmercodeiterator->fb,&charcode,err);
      if (retval < 0)
      {
        *kmercodeptr = NULL;
        return -1;
      }
      if (retval != 0)
      {
        shiftrightwithchar(kmercodeiterator->spwp,charcode);
        kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
        kmercodeiterator->currentposition++,
        *kmercodeptr = &kmercodeiterator->kmercode;
        return 0;
      }
      kmercodeiterator->inputexhausted = true;
      kmercodeiterator->totallength = kmercodeiterator->currentposition;
    } else
    {
      kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
      kmercodeiterator->hasprocessedfirst = true;
      *kmercodeptr = &kmercodeiterator->kmercode;
      return 0;
    }
  }
  if (kmercodeiterator->currentposition < kmercodeiterator->totallength +
                                          kmercodeiterator->spwp->kmersize)
  {
    shiftrightwithchar(kmercodeiterator->spwp,(GtUchar) WILDCARD);
    kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
    kmercodeiterator->currentposition++,
    *kmercodeptr = &kmercodeiterator->kmercode;
  } else
  {
    *kmercodeptr = NULL;
  }
  return 0;
}

bool gt_kmercodeiterator_inputexhausted(
                              const GtKmercodeiterator *kmercodeiterator)
{
  return kmercodeiterator->inputexhausted;
}

void gt_kmercodeiterator_delete(GtKmercodeiterator *kmercodeiterator)
{
  if (kmercodeiterator == NULL)
  {
    return;
  }
  gt_encseq_reader_delete(kmercodeiterator->esr);
  kmerstream_delete(kmercodeiterator->spwp);
  gt_sequence_buffer_delete(kmercodeiterator->fb);
  gt_free(kmercodeiterator);
}

void getencseqkmers(const GtEncseq *encseq,
                    GtReadmode readmode,
                    unsigned int kmersize,
                    void(*processkmercode)(void *,
                                           unsigned long,
                                           const GtKmercode *),
                    void *processkmercodeinfo)
{
  unsigned long currentposition = 0, totallength;
  Kmerstream *spwp;
  GtUchar charcode;
  GtEncseqReader *esr;
  unsigned int numofchars, overshoot;

  totallength = gt_encseq_total_length(encseq);
  if (totallength < (unsigned long) kmersize)
  {
    return;
  }
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  spwp = kmerstream_new(numofchars,kmersize);
  esr = gt_encseq_create_reader_with_readmode(encseq,readmode,0);
  for (currentposition = 0; currentposition < (unsigned long) kmersize;
       currentposition++)
  {
    charcode = gt_encseq_reader_next_encoded_char(esr);
    GT_CHECKENCCHAR(charcode,encseq,currentposition,readmode);
    spwp->windowwidth++;
    updatespecialpositions(spwp,charcode,false,0);
    spwp->cyclicwindow[spwp->windowwidth-1] = charcode;
  }
  kmerstream_newcode(&spwp->currentkmercode,spwp);
  processkmercode(processkmercodeinfo,0,&spwp->currentkmercode);
  for (currentposition = (unsigned long) kmersize; currentposition<totallength;
       currentposition++)
  {
    charcode = gt_encseq_reader_next_encoded_char(esr);
    GT_CHECKENCCHAR(charcode,encseq,currentposition,readmode);
    shiftrightwithchar(spwp,charcode);
    kmerstream_newcode(&spwp->currentkmercode,spwp);
    processkmercode(processkmercodeinfo,currentposition + 1 - spwp->kmersize,
                    &spwp->currentkmercode);
  }
  gt_encseq_reader_delete(esr);
  for (overshoot=0; overshoot<kmersize; overshoot++)
  {
    shiftrightwithchar(spwp,(GtUchar) WILDCARD);
    kmerstream_newcode(&spwp->currentkmercode,spwp);
    processkmercode(processkmercodeinfo,
                    overshoot + currentposition + 1 - spwp->kmersize,
                    &spwp->currentkmercode);
  }
  kmerstream_delete(spwp);
}
