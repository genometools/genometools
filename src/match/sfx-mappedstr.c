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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
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
#include "stamp.h"
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
  gt_multimappowerfree(&spwp->multimappower);
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

typedef struct
{
  int shiftright;
  const GtTwobitencoding *tbptr;
  GtTwobitencoding currentencoding;
} Singlecharacterbitstreamstate;

static void scbs_init(Singlecharacterbitstreamstate *scbs,
                      const GtTwobitencoding *twobitencoding,
                      unsigned int kmersize)
{
  scbs->tbptr = twobitencoding;
  if (kmersize == 0)
  {
    scbs->currentencoding = 0;
    scbs->shiftright = 0;
  } else
  {
    scbs->currentencoding = *(scbs->tbptr++);
    gt_assert(2U * kmersize < (unsigned int) GT_INTWORDSIZE);
    scbs->shiftright = GT_INTWORDSIZE - GT_MULT2(kmersize);
  }
}

static inline GtUchar scbs_next(Singlecharacterbitstreamstate *scbs)
{
  if (scbs->shiftright > 0)
  {
    scbs->shiftright -= 2;
  } else
  {
    scbs->currentencoding = *(scbs->tbptr++);
    scbs->shiftright = GT_INTWORDSIZE-2;
  }
  return (GtUchar) (scbs->currentencoding >> scbs->shiftright) & 3;
}

typedef struct
{
  const GtTwobitencoding *tbptr;
  GtTwobitencoding currentencoding;
  unsigned int unitoffset, shiftleft, shiftright;
  GtCodetype maskright;
} Multicharacterbitstreamstate;

static void mcbs_init(Multicharacterbitstreamstate *mcbs,
                      const GtTwobitencoding *twobitencoding,
                      unsigned int kmersize)
{
  mcbs->tbptr = twobitencoding;
  mcbs->unitoffset = 0;
  mcbs->shiftleft = 2U;
  mcbs->shiftright = (unsigned int) GT_MULT2(GT_UNITSIN2BITENC - kmersize);
  mcbs->maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;
  mcbs->currentencoding = *(mcbs->tbptr++);
}

static inline GtCodetype mcbs_next(Multicharacterbitstreamstate *mcbs,
                                   unsigned int kmersize)
{
  GtCodetype kmer;

  if (mcbs->unitoffset <= (unsigned int) GT_UNITSIN2BITENC - kmersize)
  {
    kmer = (GtCodetype) (mcbs->currentencoding >> mcbs->shiftright)
                        & mcbs->maskright;
    mcbs->shiftright-=2U;
  } else
  {
    kmer = (GtCodetype)
             ((mcbs->currentencoding << mcbs->shiftleft) |
              (*(mcbs->tbptr) >>
               (GT_MULT2(GT_UNITSIN2BITENC)-mcbs->shiftleft)))
              & mcbs->maskright;
    mcbs->shiftleft+=2;
  }
  if (mcbs->unitoffset < (unsigned int) GT_UNITSIN2BITENC-1)
  {
    mcbs->unitoffset++;
  } else
  {
    mcbs->unitoffset = 0;
    mcbs->shiftleft = 2U;
    mcbs->shiftright = (unsigned int) GT_MULT2(GT_UNITSIN2BITENC - kmersize);
    mcbs->currentencoding = *(mcbs->tbptr++);
  }
  return kmer;
}

static void showdifferentkmers(int line,unsigned long pos,GtCodetype kmer1,
                               GtCodetype kmer2)
{
  char buffer[2*GT_INTWORDSIZE+1];

  fprintf(stderr,"line %d: pos=%lu\n",line,pos);
  gt_bitsequence_tostring_units(buffer,(GtBitsequence) kmer1,2U);
  fprintf(stderr,"kmer1=%s\n",buffer);
  gt_bitsequence_tostring_units(buffer,(GtBitsequence) kmer2,2U);
  fprintf(stderr,"kmer2=%s\n",buffer);
  fprintf(stderr,"kmer1=%lu != %lu= kmer2\n",kmer1,kmer2);
}

static GtCodetype gt_kmercodeatpos(const GtTwobitencoding *twobitencoding,
                                   unsigned long pos,
                                   unsigned int kmersize)
{
  unsigned int unitoffset = (unsigned int) GT_MODBYUNITSIN2BITENC(pos);
  unsigned long unitindex = GT_DIVBYUNITSIN2BITENC(pos);
  const GtCodetype maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;

  if (unitoffset <= (unsigned int) GT_UNITSIN2BITENC - kmersize)
  {
    return (GtCodetype) (twobitencoding[unitindex]
            >> GT_MULT2(GT_UNITSIN2BITENC - kmersize - unitoffset))
           & maskright;
  } else
  {
    unsigned int shiftleft = GT_MULT2(unitoffset+kmersize-GT_UNITSIN2BITENC);
    return (GtCodetype)
           ((twobitencoding[unitindex] << shiftleft) |
            (twobitencoding[unitindex+1] >> (GT_MULT2(GT_UNITSIN2BITENC) -
                                             shiftleft)))
           & maskright;
  }
}

static GtCodetype gt_kmercodefirstpos(const GtTwobitencoding *twobitencoding,
                                      unsigned int kmersize)
{
  return (GtCodetype) (twobitencoding[0] >>
                       GT_MULT2(GT_UNITSIN2BITENC - kmersize))
                       & ((GtCodetype) (1 << GT_MULT2(kmersize))-1);
}

#define GT_SWAPBITPAIRS(L1,L2,D) ((kmer & (3UL << L1)) >> D) |\
                                 ((kmer & (3UL << L2)) << D)

static GtCodetype gt_reversekmer(GtCodetype kmer,unsigned int kmersize)
{
  switch (kmersize)
  {
    case 2:
      return GT_SWAPBITPAIRS(2,0,2);
    case 3:
      return GT_SWAPBITPAIRS(4,0,4) |
             (kmer & (3U << 2));
    case 4:
      return GT_SWAPBITPAIRS(6,0,6) |
             GT_SWAPBITPAIRS(4,2,2);
    case 5:
      return GT_SWAPBITPAIRS(8,0,8) |
             GT_SWAPBITPAIRS(6,2,4) |
             (kmer & (3U << 4));
    case 6:
      return GT_SWAPBITPAIRS(10,0,10) |
             GT_SWAPBITPAIRS(8,2,6) |
             GT_SWAPBITPAIRS(6,4,2);
    case 7:
      return GT_SWAPBITPAIRS(12,0,12) |
             GT_SWAPBITPAIRS(10,2,8) |
             GT_SWAPBITPAIRS(8,4,4) |
             (kmer & (3U << 6));
    case 8:
      return GT_SWAPBITPAIRS(14,0,14) |
             GT_SWAPBITPAIRS(12,2,10) |
             GT_SWAPBITPAIRS(10,4,6) |
             GT_SWAPBITPAIRS(8,6,2);
    case 9:
      return GT_SWAPBITPAIRS(16,0,16) |
             GT_SWAPBITPAIRS(14,2,12) |
             GT_SWAPBITPAIRS(12,4,8) |
             GT_SWAPBITPAIRS(10,6,4) |
             (kmer & (3U << 8));
    case 10:
      return GT_SWAPBITPAIRS(18,0,18) |
             GT_SWAPBITPAIRS(16,2,14) |
             GT_SWAPBITPAIRS(14,4,10) |
             GT_SWAPBITPAIRS(12,6,6) |
             GT_SWAPBITPAIRS(10,8,2);
    case 11:
      return GT_SWAPBITPAIRS(20,0,20) |
             GT_SWAPBITPAIRS(18,2,16) |
             GT_SWAPBITPAIRS(16,4,12) |
             GT_SWAPBITPAIRS(14,6,8) |
             GT_SWAPBITPAIRS(12,8,4) |
             (kmer & (3U << 10));
    case 12:
      return GT_SWAPBITPAIRS(22,0,22) |
             GT_SWAPBITPAIRS(20,2,18) |
             GT_SWAPBITPAIRS(18,4,14) |
             GT_SWAPBITPAIRS(16,6,10) |
             GT_SWAPBITPAIRS(14,8,6) |
             GT_SWAPBITPAIRS(12,10,2);
    case 13:
      return GT_SWAPBITPAIRS(24,0,24) |
             GT_SWAPBITPAIRS(22,2,20) |
             GT_SWAPBITPAIRS(20,4,16) |
             GT_SWAPBITPAIRS(18,6,12) |
             GT_SWAPBITPAIRS(16,8,8) |
             GT_SWAPBITPAIRS(14,10,4) |
             (kmer & (3U << 12));
    case 14:
      return GT_SWAPBITPAIRS(26,0,26) |
             GT_SWAPBITPAIRS(24,2,22) |
             GT_SWAPBITPAIRS(22,4,18) |
             GT_SWAPBITPAIRS(20,6,14) |
             GT_SWAPBITPAIRS(18,8,10) |
             GT_SWAPBITPAIRS(16,10,6) |
             GT_SWAPBITPAIRS(14,12,2);
    default: fprintf(stderr,"illegal kmersize=%u\n",kmersize);
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static GtCodetype gt_complementkmer(GtCodetype kmer,GtCodetype maskright)
{
  return kmer ^ maskright;
}

#ifdef SKDEBUG
static void checkallreversebitpairs(void)
{
  unsigned int kmersize, code, coderev, coderevrev, maxcode;

  for (kmersize = 2U; kmersize <= 14U; kmersize++)
  {
    maxcode = (1U << 2 * kmersize)-1;
    printf("kmsize=%u,maxcode=%u\n",kmersize,maxcode);
    for (code = 0; code <= maxcode; code++)
    {
      coderev = gt_reversekmer(code,kmersize);
      coderevrev = gt_reversekmer(coderev,kmersize);
      gt_assert(coderevrev != code);
    }
  }
}
#endif

#define READNEXTCODEANDCHECKIGNORESPECIAL(POS)\
        gt_assert(kmercodeiterator != NULL);\
        kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator);\
        gt_assert(kmercodeptr != NULL);\
        if (!kmercodeptr->definedspecialposition && kmer != kmercodeptr->code)\
        {\
          showdifferentkmers(__LINE__,POS,kmer,kmercodeptr->code);\
          exit(EXIT_FAILURE);\
        }

#define READNEXTCODEANDCHECKNOSPECIAL(CODE,POS)\
        gt_assert(kmercodeiterator != NULL);\
        kmercodeptr = gt_kmercodeiterator_encseq_nonspecial_next(\
                                   kmercodeiterator);\
        gt_assert(kmercodeptr != NULL && !kmercodeptr->definedspecialposition);\
        if ((CODE) != kmercodeptr->code)\
        {\
          showdifferentkmers(__LINE__,pos,CODE,kmercodeptr->code);\
          exit(EXIT_FAILURE);\
        }

#define UPDATEKMER(KMER,CC)\
        KMER <<= 2;\
        KMER |= cc;\
        KMER &= maskright

static GtCodetype getencseqkmers_nospecialtwobitencoding(
                                            const GtEncseq *encseq,
                                            GtReadmode readmode,
                                            unsigned int kmersize,
                                            void(*processkmercode)(void *,
                                                           unsigned long,
                                                           GtCodetype),
                                            void *processkmercodeinfo,
                                            unsigned long startpos,
                                            unsigned long endpos,
                                            GtKmercodeiterator
                                              *kmercodeiterator)
{
  unsigned long pos, unitindex, totallength;
  GtKmercode kmer;
  GtUchar cc;
  const GtTwobitencoding *twobitencoding;
  GtTwobitencoding currentencoding;
  unsigned int shiftright;
  GtCodetype complcode;
  const GtCodetype maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;
  const GtKmercode *kmercodeptr;

  twobitencoding = gt_encseq_twobitencoding_export(encseq);
  totallength = gt_encseq_total_length(encseq);
  kmer.definedspecialposition = false;
  kmer.specialposition = 0;
  if (GT_ISDIRREVERSE(readmode))
  {
    gt_assert(endpos >= (unsigned long) kmersize);
    pos = endpos - (unsigned long) kmersize;
    unitindex = (pos > 0) ? GT_DIVBYUNITSIN2BITENC(pos-1) : 0;
    if (kmersize > 1U)
    {
      kmer.code = gt_reversekmer(gt_kmercodeatpos(twobitencoding,pos,kmersize),
                                 kmersize);
    }
    if (readmode == GT_READMODE_REVCOMPL)
    {
      complcode = gt_complementkmer(kmer.code,maskright);
      READNEXTCODEANDCHECKNOSPECIAL(complcode,pos);
      processkmercode(processkmercodeinfo,pos,complcode);
    } else
    {
      READNEXTCODEANDCHECKNOSPECIAL(kmer.code,pos);
      processkmercode(processkmercodeinfo,pos,kmer.code);
    }
  } else
  {
    pos = startpos;
    unitindex = GT_DIVBYUNITSIN2BITENC(startpos+kmersize);
    kmer.code = gt_kmercodeatpos(twobitencoding,pos,kmersize);
    if (readmode == GT_READMODE_COMPL)
    {
      complcode = gt_complementkmer(kmer.code,maskright);
      READNEXTCODEANDCHECKNOSPECIAL(complcode,pos);
      processkmercode(processkmercodeinfo,pos,complcode);
    } else
    {
      READNEXTCODEANDCHECKNOSPECIAL(kmer.code,pos);
      processkmercode(processkmercodeinfo,pos,kmer.code);
    }
  }
  currentencoding = twobitencoding[unitindex];
  if (GT_ISDIRREVERSE(readmode))
  {
    shiftright = (unsigned int)
                 GT_MULT2(GT_UNITSIN2BITENC - 1 -
                          GT_MODBYUNITSIN2BITENC(pos-1));
    while (pos > startpos)
    {
      pos--;
      cc = (GtUchar) (currentencoding >> shiftright) & 3;
      UPDATEKMER(kmer.code,cc);
      if (readmode == GT_READMODE_REVCOMPL)
      {
        complcode = gt_complementkmer(kmer.code,maskright);
        READNEXTCODEANDCHECKNOSPECIAL(complcode,pos);
        processkmercode(processkmercodeinfo,pos,complcode);
      } else
      {
        READNEXTCODEANDCHECKNOSPECIAL(kmer.code,pos);
        processkmercode(processkmercodeinfo,pos,kmer.code);
      }
      if (shiftright < (unsigned int) (GT_INTWORDSIZE-2))
      {
        shiftright += 2;
      } else
      {
        if (unitindex > 0)
        {
          currentencoding = twobitencoding[--unitindex];
        }
        shiftright = 0;
      }
    }
  } else
  {
    unsigned long maxunitindex = gt_unitsoftwobitencoding(totallength) - 1;

    shiftright = (unsigned int)
                 GT_MULT2(GT_UNITSIN2BITENC - 1 -
                          GT_MODBYUNITSIN2BITENC(startpos+kmersize));
    while (pos < endpos - (unsigned long) kmersize)
    {
      pos++;
      cc = (GtUchar) (currentencoding >> shiftright) & 3;
      UPDATEKMER(kmer.code,cc);
      if (readmode == GT_READMODE_COMPL)
      {
        complcode = gt_complementkmer(kmer.code,maskright);
        READNEXTCODEANDCHECKNOSPECIAL(complcode,pos);
        processkmercode(processkmercodeinfo,pos,complcode);
      } else
      {
        READNEXTCODEANDCHECKNOSPECIAL(kmer.code,pos);
        processkmercode(processkmercodeinfo,pos,kmer.code);
      }
      if (shiftright > 0)
      {
        shiftright -= 2;
      } else
      {
        if (unitindex < maxunitindex-1)
        {
          currentencoding = twobitencoding[++unitindex];
        }
        shiftright = (unsigned int) (GT_INTWORDSIZE-2);
      }
    }
  }
  return kmer.code;
}

static unsigned long getencseqkmers_rangetwobitencoding(const GtEncseq *encseq,
                                               GtReadmode readmode,
                                               unsigned int kmersize,
                                               void(*processkmercode)(void *,
                                                           unsigned long,
                                                           GtCodetype),
                                               void *processkmercodeinfo,
                                               Codeatposition *codelist,
                                               unsigned long nextspecialcode,
                                               unsigned long startpos,
                                               unsigned long endpos,
                                               GtKmercodeiterator
                                                 *kmercodeiterator)
{
  GtCodetype lastcode;
  const GtCodetype maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;
  GtCodetype newcode;
  unsigned long totallength = gt_encseq_total_length(encseq);

  if (endpos - startpos >= (unsigned long) kmersize)
  {
    gt_assert(endpos > 0);
    lastcode = getencseqkmers_nospecialtwobitencoding(encseq,
                                                      readmode,
                                                      kmersize,
                                                      processkmercode,
                                                      processkmercodeinfo,
                                                      startpos,
                                                      endpos,
                                                      kmercodeiterator);
    if (GT_ISDIRCOMPLEMENT(readmode))
    {
      lastcode = gt_complementkmer(lastcode,maskright);
    }
    newcode = ((lastcode << 2) | 3UL) & maskright;
    if (codelist != NULL)
    {
      codelist[nextspecialcode].maxprefixindex = kmersize - 1;
      codelist[nextspecialcode].code = (unsigned int) newcode;
      codelist[nextspecialcode].position
        = GT_ISDIRREVERSE(readmode) ? (totallength - startpos) : endpos;
      /*
      printf("fast1: store(code=%u,maxprefixindex=%u,pos=%lu)\n",
               codelist[nextspecialcode].code,
               codelist[nextspecialcode].maxprefixindex,
               codelist[nextspecialcode].position);
      */
    }
    return nextspecialcode + 1;
  }
  if (startpos < endpos)
  {
    unsigned int fillpos;
    const GtTwobitencoding *twobitencoding;

    gt_assert((unsigned long) kmersize > endpos - startpos);
    fillpos = (unsigned int) (kmersize - (endpos - startpos));
    twobitencoding = gt_encseq_twobitencoding_export(encseq);
    lastcode = gt_kmercodeatpos(twobitencoding,startpos,
                                (unsigned int) (endpos - startpos));
    if (GT_ISDIRREVERSE(readmode) && (unsigned int) (endpos - startpos) > 1U)
    {
      lastcode = gt_reversekmer(lastcode,(unsigned int) (endpos - startpos));
    }
    if (GT_ISDIRCOMPLEMENT(readmode))
    {
      lastcode = gt_complementkmer(lastcode,maskright);
    }
    newcode
      = ((lastcode << GT_MULT2(fillpos)) | ((1UL << GT_MULT2(fillpos)) - 1))
         & maskright;
    if (codelist != NULL)
    {
      codelist[nextspecialcode].maxprefixindex
        = (unsigned int) (endpos - startpos);
      codelist[nextspecialcode].code = (unsigned int) newcode;
      codelist[nextspecialcode].position
        = GT_ISDIRREVERSE(readmode) ? (totallength - startpos) : endpos;
      /*
      printf("fast2: store(code=%u,maxprefixindex=%u,pos=%lu)\n",
               codelist[nextspecialcode].code,
               codelist[nextspecialcode].maxprefixindex,
               codelist[nextspecialcode].position);
      */
    }
    return nextspecialcode + 1;
  }
  return nextspecialcode;
}

unsigned long getencseqkmers_twobitencoding(
                                   const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   unsigned int kmersize,
                                   void(*processkmercode)(void *,
                                                          unsigned long,
                                                          GtCodetype),
                                   void *processkmercodeinfo,
                                   Codeatposition *codelist)
{
  unsigned long nextspecialcode = 0, laststart = 0, lastend, totallength;
  GtKmercodeiterator *kmercodeiterator
    = gt_kmercodeiterator_encseq_new(encseq,readmode,kmersize,0);

  lastend = totallength = gt_encseq_total_length(encseq);
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;

    if (GT_ISDIRREVERSE(readmode))
    {
      sri = gt_specialrangeiterator_new(encseq,false);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        gt_assert(range.end <= lastend);
        nextspecialcode
          = getencseqkmers_rangetwobitencoding(encseq,
                                               readmode,
                                               kmersize,
                                               processkmercode,
                                               processkmercodeinfo,
                                               codelist,
                                               nextspecialcode,
                                               range.end,
                                               lastend,
                                               kmercodeiterator);
        lastend = range.start;
      }
    } else
    {
      sri = gt_specialrangeiterator_new(encseq,true);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        gt_assert(range.start >= laststart);
        nextspecialcode
          = getencseqkmers_rangetwobitencoding(encseq,
                                               readmode,
                                               kmersize,
                                               processkmercode,
                                               processkmercodeinfo,
                                               codelist,
                                               nextspecialcode,
                                               laststart,
                                               range.start,
                                               kmercodeiterator);
        laststart = range.end;
      }
    }
    gt_assert(totallength >= laststart);
    gt_specialrangeiterator_delete(sri);
  }
  if (GT_ISDIRREVERSE(readmode))
  {
    nextspecialcode = getencseqkmers_rangetwobitencoding(encseq,
                                                         readmode,
                                                         kmersize,
                                                         processkmercode,
                                                         processkmercodeinfo,
                                                         codelist,
                                                         nextspecialcode,
                                                         0,
                                                         lastend,
                                                         kmercodeiterator);
  } else
  {
    nextspecialcode = getencseqkmers_rangetwobitencoding(encseq,
                                                         readmode,
                                                         kmersize,
                                                         processkmercode,
                                                         processkmercodeinfo,
                                                         codelist,
                                                         nextspecialcode,
                                                         laststart,
                                                         totallength,
                                                         kmercodeiterator);
  }
  gt_kmercodeiterator_delete(kmercodeiterator);
  return nextspecialcode;
}

static void swallowkmercode(GT_UNUSED void *processinfo,
                            GT_UNUSED unsigned long pos,
                            GT_UNUSED GtCodetype code)
{
  return;
}

static void gt_encseq_faststream_kmers(const GtEncseq *encseq,
                                       Bitstreamreadmode bsrsmode,
                                       unsigned int kmersize)
{
  unsigned long totallength, pos;
  GtCodetype kmer;
  GtKmercodeiterator *kmercodeiterator = NULL;
  const GtKmercode *kmercodeptr;
  const GtTwobitencoding *twobitencoding;
  Multicharacterbitstreamstate mcbs;

  gt_assert(kmersize < (unsigned int) GT_UNITSIN2BITENC);
  totallength = gt_encseq_total_length(encseq);
  if (totallength < (unsigned long) kmersize)
  {
    return;
  }
  twobitencoding = gt_encseq_twobitencoding_export(encseq);
  if (bsrsmode == BSRS_reader_multi ||
      bsrsmode == BSRS_stream_reader_multi ||
      bsrsmode == BSRS_stream_reader_multi2)
  {
    kmercodeiterator = gt_kmercodeiterator_encseq_new(encseq,
                                                      GT_READMODE_FORWARD,
                                                      kmersize,0);
  }
  switch (bsrsmode)
  {
    case BSRS_stream_multi:
      {
        uint64_t kmersum = 0;
        GtCodetype kmer2;
        mcbs_init(&mcbs,twobitencoding,kmersize);
        for (pos = 0; pos <= totallength - (unsigned long) kmersize; pos++)
        {
          kmer = mcbs_next(&mcbs,kmersize);
          kmer2 = gt_kmercodeatpos(twobitencoding, pos, kmersize);
          if (kmer != kmer2)
          {
            showdifferentkmers(__LINE__,pos,kmer,kmer2);
          }
          gt_assert(kmer == kmer2);
          kmersum += (uint64_t) kmer;
        }
        printf("kmersum=" Formatuint64_t "\n",PRINTuint64_tcast(kmersum));
      }
      break;
    case BSRS_reader_multi:
      {
        uint64_t kmersum = 0;

        for (pos = 0; pos <= totallength - (unsigned long) kmersize; pos++)
        {
          kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator);
          gt_assert(kmercodeptr != NULL);
          kmersum += (uint64_t) kmercodeptr->code;
        }
        printf("kmersum=" Formatuint64_t "\n",PRINTuint64_tcast(kmersum));
        break;
      }
    case BSRS_stream_reader_multi:
      mcbs_init(&mcbs,twobitencoding,kmersize);
      for (pos = 0; pos <= totallength - (unsigned long) kmersize; pos++)
      {
        kmer = mcbs_next(&mcbs,kmersize);
        READNEXTCODEANDCHECKIGNORESPECIAL(pos);
      }
      break;
    case BSRS_stream_reader_multi2:
      {
        Singlecharacterbitstreamstate scbs;
        GtUchar cc;
        const GtCodetype maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;
        uint64_t kmersum = 0;

        kmer = gt_kmercodefirstpos(twobitencoding,kmersize);
        kmersum += (uint64_t) kmer;
        scbs_init(&scbs,twobitencoding,kmersize);
        READNEXTCODEANDCHECKIGNORESPECIAL(0);
        for (pos = 1UL; pos <= totallength - (unsigned long) kmersize; pos++)
        {
          cc = scbs_next(&scbs);
          UPDATEKMER(kmer,cc);
          READNEXTCODEANDCHECKIGNORESPECIAL(pos);
          kmersum += (uint64_t) kmer;
        }
        printf("kmersum=" Formatuint64_t "\n",PRINTuint64_tcast(kmersum));
        break;
      }
    case BSRS_stream_reader_multi3:
      printf("getencseqkmers_twobitencoding(kmersize=%u,forward)\n",kmersize);
      (void) getencseqkmers_twobitencoding(encseq,
                                           GT_READMODE_FORWARD,
                                           kmersize,
                                           swallowkmercode,
                                           NULL,NULL);
      printf("getencseqkmers_twobitencoding(kmersize=%u,reverse)\n",kmersize);
      (void) getencseqkmers_twobitencoding(encseq,
                                           GT_READMODE_REVERSE,
                                           kmersize,
                                           swallowkmercode,
                                           NULL,NULL);
      printf("getencseqkmers_twobitencoding(kmersize=%u,compl)\n",kmersize);
      (void) getencseqkmers_twobitencoding(encseq,
                                           GT_READMODE_COMPL,
                                           kmersize,
                                           swallowkmercode,
                                           NULL,NULL);
      printf("getencseqkmers_twobitencoding(kmersize=%u,revcompl)\n",kmersize);
      (void) getencseqkmers_twobitencoding(encseq,
                                           GT_READMODE_REVCOMPL,
                                           kmersize,
                                           swallowkmercode,
                                           NULL,NULL);
      break;
    default:
      break;
  }
  gt_kmercodeiterator_delete(kmercodeiterator);
}

void gt_encseq_faststream(const GtEncseq *encseq,
                          Bitstreamreadmode bsrsmode,
                          unsigned int multiarg)
{
  const GtTwobitencoding *twobitencoding;

  twobitencoding = gt_encseq_twobitencoding_export(encseq);
  if (twobitencoding != NULL)
  {
    unsigned long idx, totallength, pos;
    uint64_t pairbitsum = 0, pairbitsumBF;
    GtUchar cc, ccesr;
    GtEncseqReader *esr = NULL;
    Singlecharacterbitstreamstate scbs;

    scbs_init(&scbs,twobitencoding,0);
    if (bsrsmode == BSRS_reader_single ||
        bsrsmode == BSRS_stream_reader_single)
    {
      esr = gt_encseq_create_reader_with_readmode(encseq,
                                                  GT_READMODE_FORWARD,
                                                  0);
    }
    totallength = gt_encseq_total_length(encseq);
    switch (bsrsmode)
    {
      case BSRS_stream_words:
        for (idx = 0; idx < gt_unitsoftwobitencoding(totallength); idx++)
        {
          pairbitsum += twobitencoding[idx];
        }
        break;
      case BSRS_stream_single:
        for (pos = 0; pos < totallength; pos++)
        {
          cc = scbs_next(&scbs);
          pairbitsum += (uint64_t) cc;
        }
        pairbitsumBF = gt_encseq_pairbitsum(encseq);
        if (pairbitsum != pairbitsumBF)
        {
          fprintf(stderr,"pairbitsum=" Formatuint64_t "!=" Formatuint64_t
                         "=pairbitsumBF\n",
                         PRINTuint64_tcast(pairbitsum),
                         PRINTuint64_tcast(pairbitsumBF));
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        break;
      case BSRS_reader_single:
        for (pos = 0; pos < totallength; pos++)
        {
          ccesr = gt_encseq_reader_next_encoded_char(esr);
          pairbitsum += (uint64_t) ccesr;
        }
        break;
      case BSRS_stream_reader_single:
        for (pos = 0; pos < totallength; pos++)
        {
          cc = scbs_next(&scbs);
          pairbitsum += (uint64_t) cc;
          ccesr = gt_encseq_reader_next_encoded_char(esr);
          pairbitsum += (uint64_t) ccesr;
          gt_assert(cc == ccesr || ISSPECIAL(ccesr));
        }
        break;
      case BSRS_stream_multi:
      case BSRS_reader_multi:
      case BSRS_stream_reader_multi:
      case BSRS_stream_reader_multi2:
      case BSRS_stream_reader_multi3:
        gt_encseq_faststream_kmers(encseq,bsrsmode,multiarg);
        break;
    }
    printf("pairbitsum=" Formatuint64_t "\n",PRINTuint64_tcast(pairbitsum));
    gt_encseq_reader_delete(esr);
  }
}
