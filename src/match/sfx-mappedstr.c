/*
  Copyright (c) 2007-2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "core/alphabet.h"
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/encseq.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/types_api.h"
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

static unsigned int determinefirstspecialposition(unsigned int windowwidth,
                                                  unsigned int kmersize,
                                                  const GtUchar *cyclicwindow,
                                                  unsigned int firstindex)
{
  unsigned int i;

  for (i=0; i < windowwidth; i++)
  {
    if (ISSPECIAL(cyclicwindow[(firstindex+i) % kmersize]))
    {
      return i;
    }
  }
  return kmersize;
}
#endif

/* Each special character occurring in the current window defines
   a value of type Specialcontext. The context of a special symbol N is
   the maximal sequence in the current window to the left of N consisting of
   non-special symbols only. So the context either starts with the
   first symbol in the window or right after a special symbol. We store the
   length of the context and the integer code for the context such that
   the first symbol of the context is weighted by $\sigma^{q-1}$. */

typedef struct
{
  unsigned int lengthofleftcontext;
  GtCodetype codeofleftcontext;
} GtSpecialcontext;

/* We store the special context in a queue, such that the first
   special symbol in the window is at the head of the queue and the
   last special symbol in the window is at the tail of the queue. */

typedef struct
{
  GtSpecialcontext *queuespace,  /* the space to store the queue elements */
                   *enqueueptr,  /* entry into which element is to be enqued */
                   *dequeueptr;  /* last element of queue */
  unsigned int queuesize,        /* size of the queue */
               noofelements;     /* no ofelements between enqueueptr+1 and
                                    dequeptr */
} GtSpecialqueue;

typedef struct
{
  GtSpecialqueue specialqueue;
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
} GtKmerstream;

#ifdef SKDEBUG
static void special_queue_verify(const GtKmerstream *kmerstream)
{
  unsigned int i, contextsize = 0, verified = 0;
  GtSpecialcontext *queueptr = kmerstream->specialqueue.dequeueptr;

  gt_assert(kmerstream->kmersize == kmerstream->windowwidth);
  for (i=0; i < kmerstream->kmersize; i++)
  {
    GtUchar cc = kmerstream->cyclicwindow[(kmerstream->firstindex+i) %
                                           kmerstream->kmersize];
    if (ISSPECIAL(cc))
    {
      gt_assert (queueptr->lengthofleftcontext == contextsize);
      if (queueptr > kmerstream->specialqueue.queuespace)
      {
        queueptr--;
      } else
      {
        queueptr = kmerstream->specialqueue.queuespace +
                   kmerstream->specialqueue.queuesize - 1;
      }
      verified++;
      contextsize = 0;
    } else
    {
      contextsize++;
    }
  }
  gt_assert(verified == kmerstream->specialqueue.noofelements);
}
#endif

static void special_queue_reset(GtSpecialqueue *specialqueue)
{
  specialqueue->noofelements = 0;
  specialqueue->dequeueptr = specialqueue->enqueueptr
                           = specialqueue->queuespace +
                             specialqueue->queuesize - 1;
}

static void special_queue_init(GtSpecialqueue *specialqueue,
                               unsigned int queuesize)
{
  specialqueue->queuesize = queuesize;
  specialqueue->queuespace
    = gt_malloc(queuesize * sizeof (*specialqueue->queuespace));
  special_queue_reset(specialqueue);
}

static bool special_queue_is_empty(const GtSpecialqueue *specialqueue)
{
  return (specialqueue->noofelements == 0) ? true : false;
}

static GtSpecialcontext *special_queue_head_get(const GtSpecialqueue
                                                 *specialqueue)
{
  return specialqueue->dequeueptr;
}

static void special_queueu_head_delete(GtSpecialqueue *specialqueue)
{
  specialqueue->noofelements--;
  if (specialqueue->dequeueptr > specialqueue->queuespace)
  {
    specialqueue->dequeueptr--;
  } else
  {
    specialqueue->dequeueptr
      = specialqueue->queuespace + specialqueue->queuesize - 1;
  }
}

static void special_queue_enqueue(GtSpecialqueue *specialqueue,
                                  unsigned int lengthofleftcontext,
                                  GtCodetype codeofleftcontext)
{
  specialqueue->noofelements++;
  specialqueue->enqueueptr->codeofleftcontext = codeofleftcontext;
  specialqueue->enqueueptr->lengthofleftcontext = lengthofleftcontext;
  if (specialqueue->enqueueptr > specialqueue->queuespace)
  {
    specialqueue->enqueueptr--;
  } else
  {
    specialqueue->enqueueptr
      = specialqueue->queuespace + specialqueue->queuesize - 1;
  }
}

static void special_queue_delete(GtSpecialqueue *specialqueue)
{
  gt_free(specialqueue->queuespace);
}

static void kmerstream_init(GtKmerstream *spwp,
                            unsigned int numofchars,
                            unsigned int kmersize)
{
  spwp->lengthwithoutspecial = 0;
  spwp->codewithoutspecial = 0;
  spwp->kmersize = kmersize;
  spwp->numofchars = numofchars;
  spwp->windowwidth = 0;
  spwp->firstindex = 0;
}

static GtKmerstream *kmerstream_new(unsigned int numofchars,
                                    unsigned int kmersize)
{
  GtKmerstream *spwp;

  gt_assert(kmersize <= (unsigned int) MAXPREFIXLENGTH);
  spwp = gt_malloc(sizeof (*spwp));
  spwp->multimappower = gt_initmultimappower(numofchars,kmersize);
  kmerstream_init(spwp, numofchars, kmersize);
  special_queue_init(&spwp->specialqueue,kmersize);
  spwp->filltable = gt_filllargestchartable(numofchars,kmersize);
  return spwp;
}

static void kmerstream_reset(GtKmerstream *spwp)
{
  kmerstream_init(spwp, spwp->numofchars, spwp->kmersize);
  special_queue_reset(&(spwp->specialqueue));
}

static void kmerstream_updatespecialpositions(GtKmerstream *spwp,
                                              GtUchar charcode,
                                              bool doshift,
                                              GtUchar leftchar)
{
  if (doshift && !special_queue_is_empty(&spwp->specialqueue))
  {
    GtSpecialcontext *head = special_queue_head_get(&spwp->specialqueue);

    if (head->lengthofleftcontext > 0)
    {
      SUBTRACTLCHARANDSHIFT(head->codeofleftcontext,leftchar,spwp->numofchars,
                            spwp->multimappower[0]);
      head->lengthofleftcontext--;
    } else
    {
      special_queueu_head_delete(&spwp->specialqueue);
    }
  }
  if (ISNOTSPECIAL(charcode))
  {
    if (spwp->lengthwithoutspecial == spwp->kmersize)
    {
      SUBTRACTLCHARSHIFTADDNEXT(spwp->codewithoutspecial,
                                leftchar,
                                spwp->numofchars,
                                spwp->multimappower[0],
                                charcode);
    } else
    {
      spwp->codewithoutspecial
        += spwp->multimappower[spwp->lengthwithoutspecial][charcode];
      spwp->lengthwithoutspecial++;
    }
  } else
  {
    /* only here we add some element to the queue */
    unsigned int newelem_lengthofleftcontext
      = special_queue_is_empty(&spwp->specialqueue)
          ? (spwp->windowwidth - 1U)
          : spwp->lengthwithoutspecial;
    if (spwp->lengthwithoutspecial == spwp->kmersize)
    {
      SUBTRACTLCHARANDSHIFT(spwp->codewithoutspecial,
                            leftchar,
                            spwp->numofchars,
                            spwp->multimappower[0]);
    }
    /* only here we add some element to the queue */
    gt_assert(newelem_lengthofleftcontext < spwp->kmersize);
    special_queue_enqueue(&spwp->specialqueue,
                          newelem_lengthofleftcontext,
                          spwp->codewithoutspecial);
    spwp->lengthwithoutspecial = 0;
    spwp->codewithoutspecial = 0;
  }
}

static void kmerstream_newcode(GtKmercode *kmercode, GtKmerstream *spwp)
{
  {
    if (special_queue_is_empty(&spwp->specialqueue))
    {
      kmercode->definedspecialposition = false;
      kmercode->specialposition = 0;
      kmercode->code = spwp->codewithoutspecial;
    } else
    {
      GtSpecialcontext *head = special_queue_head_get(&spwp->specialqueue);
      gt_assert(head->lengthofleftcontext < spwp->kmersize);
      kmercode->code = head->codeofleftcontext +
                       spwp->filltable[head->lengthofleftcontext];
      kmercode->definedspecialposition = true;
      kmercode->specialposition = head->lengthofleftcontext;
    }
  }
}

static void kmerstream_shiftrightwithchar(GtKmerstream *spwp, GtUchar charcode)
{
  gt_assert(spwp->windowwidth == spwp->kmersize);
  kmerstream_updatespecialpositions(spwp,charcode,true,
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

static void kmerstream_delete(GtKmerstream *spwp)
{
  if (spwp != NULL) {
    gt_free(spwp->filltable);
    gt_multimappower_delete(spwp->multimappower);
    special_queue_delete(&spwp->specialqueue);
    gt_free(spwp);
  }
}

struct GtKmercodeiterator
{
  GtUword totallength, startpos, currentposition;
  bool hasprocessedfirst, inputexhausted;
  const GtEncseq *encseq;
  GtEncseqReader *esr;
  GtReadmode readmode;
  GtKmerstream *spwp;
  GtKmercode kmercode;
  GtSequenceBuffer *fb; /* only for generating from file */
};

void gt_kmercodeiterator_reset(GtKmercodeiterator *kmercodeiterator,
                               GtReadmode readmode,
                               GtUword startpos)
{
  GtUchar charcode;
  const GtEncseq *encseq = kmercodeiterator->encseq;
  GtUword kmersize = (GtUword) kmercodeiterator->spwp->kmersize;

  gt_assert(!GT_ISDIRREVERSE(readmode) || startpos == 0);
  kmercodeiterator->totallength = gt_encseq_total_length(encseq);
  kmercodeiterator->startpos = startpos;
  gt_assert(startpos < kmercodeiterator->totallength);
  kmercodeiterator->fb = NULL;
  if (kmercodeiterator->totallength - startpos < kmersize)
  {
    kmercodeiterator->inputexhausted = true;
    gt_encseq_reader_delete(kmercodeiterator->esr);
    kmercodeiterator->esr = NULL;
    kmerstream_delete(kmercodeiterator->spwp);
    kmercodeiterator->spwp = NULL;
  } else
  {
    kmercodeiterator->inputexhausted = false;
    kmercodeiterator->readmode = readmode;
    gt_encseq_reader_reinit_with_readmode(kmercodeiterator->esr,
                                          encseq,
                                          readmode,
                                          startpos);
    kmerstream_reset(kmercodeiterator->spwp);
    kmercodeiterator->hasprocessedfirst = false;
    for (kmercodeiterator->currentposition = startpos;
         kmercodeiterator->currentposition < startpos+(GtUword) kmersize;
         kmercodeiterator->currentposition++)
    {
      charcode = gt_encseq_reader_next_encoded_char(kmercodeiterator->esr);
      kmercodeiterator->spwp->windowwidth++;
      kmerstream_updatespecialpositions(kmercodeiterator->spwp,charcode,
                                        false,0);
      kmercodeiterator->spwp->cyclicwindow[kmercodeiterator->
                                           spwp->windowwidth-1] = charcode;
    }
  }
}

/*@notnull@*/ GtKmercodeiterator *gt_kmercodeiterator_encseq_new(
                                            const GtEncseq *encseq,
                                            GtReadmode readmode,
                                            unsigned int kmersize,
                                            GtUword startpos)
{
  GtKmercodeiterator *kmercodeiterator;
  unsigned int numofchars;
  GtUchar charcode;

  gt_assert(!GT_ISDIRREVERSE(readmode) || startpos == 0);
  kmercodeiterator = gt_malloc(sizeof (*kmercodeiterator));
  kmercodeiterator->totallength = gt_encseq_total_length(encseq);
  kmercodeiterator->startpos = startpos;
  gt_assert(startpos < kmercodeiterator->totallength);
  kmercodeiterator->fb = NULL;
  kmercodeiterator->encseq = encseq;
  if (kmercodeiterator->totallength - startpos < (GtUword) kmersize)
  {
    kmercodeiterator->inputexhausted = true;
    kmercodeiterator->esr = NULL;
    kmercodeiterator->spwp = NULL;
  } else
  {
    kmercodeiterator->inputexhausted = false;
    kmercodeiterator->readmode = readmode;
    kmercodeiterator->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                                  readmode,
                                                                  startpos);
    numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
    kmercodeiterator->spwp = kmerstream_new(numofchars,kmersize);
    kmercodeiterator->hasprocessedfirst = false;
    for (kmercodeiterator->currentposition = startpos;
         kmercodeiterator->currentposition < startpos+(GtUword) kmersize;
         kmercodeiterator->currentposition++)
    {
      charcode = gt_encseq_reader_next_encoded_char(kmercodeiterator->esr);
      kmercodeiterator->spwp->windowwidth++;
      kmerstream_updatespecialpositions(kmercodeiterator->spwp,charcode,
                                        false,0);
      kmercodeiterator->spwp->cyclicwindow[kmercodeiterator->
                                           spwp->windowwidth-1] = charcode;
    }
  }
  return kmercodeiterator;
}

/* TODO: remove this. Outside world should handle this on its own. */
GtUword
gt_kmercodeiterator_encseq_get_currentpos(GtKmercodeiterator *kmercodeiterator)
{
  return kmercodeiterator->currentposition;
}

void
gt_kmercodeiterator_encseq_setexhausted(GtKmercodeiterator *kmercodeiterator,
                                        bool exhausted)
{
  kmercodeiterator->inputexhausted = exhausted;
}

const GtKmercode *gt_kmercodeiterator_encseq_next(
                             GtKmercodeiterator *kmercodeiterator)
{
  if (!kmercodeiterator->hasprocessedfirst)
  {
    gt_assert(kmercodeiterator->currentposition
              == kmercodeiterator->startpos + kmercodeiterator->spwp->kmersize);
    kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
    kmercodeiterator->hasprocessedfirst = true;
    return &kmercodeiterator->kmercode;
  }
  if (kmercodeiterator->currentposition < kmercodeiterator->totallength)
  {
    GtUchar charcode
      = gt_encseq_reader_next_encoded_char(kmercodeiterator->esr);
    kmerstream_shiftrightwithchar(kmercodeiterator->spwp,charcode);
    kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
    kmercodeiterator->currentposition++;
    return &kmercodeiterator->kmercode;
  }
  if (kmercodeiterator->currentposition < kmercodeiterator->totallength +
                                          kmercodeiterator->spwp->kmersize)
  {
    kmerstream_shiftrightwithchar(kmercodeiterator->spwp,(GtUchar) WILDCARD);
    kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
    kmercodeiterator->currentposition++;
    return &kmercodeiterator->kmercode;
  }
  return NULL;
}

const GtKmercode *gt_kmercodeiterator_encseq_nonspecial_next(
                             GtKmercodeiterator *kmercodeiterator)
{
  while (true) {
    if (!kmercodeiterator->hasprocessedfirst) {
      gt_assert(kmercodeiterator->currentposition
                == kmercodeiterator->startpos +
                   kmercodeiterator->spwp->kmersize);
      kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
      kmercodeiterator->hasprocessedfirst = true;
      if (!kmercodeiterator->kmercode.definedspecialposition)
        return &kmercodeiterator->kmercode;
    }
    else {
      if (kmercodeiterator->currentposition < kmercodeiterator->totallength) {
        GtUchar charcode
          = gt_encseq_reader_next_encoded_char(kmercodeiterator->esr);
        kmerstream_shiftrightwithchar(kmercodeiterator->spwp,charcode);
        kmerstream_newcode(&kmercodeiterator->kmercode, kmercodeiterator->spwp);
        kmercodeiterator->currentposition++;
        if (!kmercodeiterator->kmercode.definedspecialposition)
          return &kmercodeiterator->kmercode;
      }
      else
        break;
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
         kmercodeiterator->currentposition < (GtUword) kmersize;
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
      kmerstream_updatespecialpositions(kmercodeiterator->spwp,charcode,
                                        false,0);
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
        kmerstream_shiftrightwithchar(kmercodeiterator->spwp,charcode);
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
    kmerstream_shiftrightwithchar(kmercodeiterator->spwp,(GtUchar) WILDCARD);
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
  if (kmercodeiterator != NULL)
  {
    gt_encseq_reader_delete(kmercodeiterator->esr);
    kmerstream_delete(kmercodeiterator->spwp);
    gt_sequence_buffer_delete(kmercodeiterator->fb);
    gt_free(kmercodeiterator);
  }
}

void getencseqkmers(const GtEncseq *encseq,
                    GtReadmode readmode,
                    unsigned int kmersize,
                    void(*processkmercode)(void *,
                                           GtUword,
                                           const GtKmercode *),
                    void *processkmercodeinfo)
{
  GtUword currentposition = 0, totallength;
  GtKmerstream *spwp;
  GtUchar charcode;
  GtEncseqReader *esr;
  unsigned int numofchars, overshoot;

  totallength = gt_encseq_total_length(encseq);
  if (totallength < (GtUword) kmersize)
  {
    return;
  }
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  spwp = kmerstream_new(numofchars,kmersize);
  esr = gt_encseq_create_reader_with_readmode(encseq,readmode,0);
  for (currentposition = 0; currentposition < (GtUword) kmersize;
       currentposition++)
  {
    charcode = gt_encseq_reader_next_encoded_char(esr);
    spwp->windowwidth++;
    kmerstream_updatespecialpositions(spwp,charcode,false,0);
    spwp->cyclicwindow[spwp->windowwidth-1] = charcode;
  }
  kmerstream_newcode(&spwp->currentkmercode,spwp);
  processkmercode(processkmercodeinfo,0,&spwp->currentkmercode);
  for (currentposition = (GtUword) kmersize; currentposition<totallength;
       currentposition++)
  {
    charcode = gt_encseq_reader_next_encoded_char(esr);
    kmerstream_shiftrightwithchar(spwp,charcode);
    kmerstream_newcode(&spwp->currentkmercode,spwp);
    processkmercode(processkmercodeinfo,currentposition + 1 - spwp->kmersize,
                    &spwp->currentkmercode);
  }
  gt_encseq_reader_delete(esr);
  for (overshoot=0; overshoot<kmersize; overshoot++)
  {
    kmerstream_shiftrightwithchar(spwp,(GtUchar) WILDCARD);
    kmerstream_newcode(&spwp->currentkmercode,spwp);
    processkmercode(processkmercodeinfo,
                    overshoot + currentposition + 1 - spwp->kmersize,
                    &spwp->currentkmercode);
  }
  kmerstream_delete(spwp);
}
