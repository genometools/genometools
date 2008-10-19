/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "encseq-def.h"
#include "iter-window.h"

struct Windowiterator
{
  Uchar *buffer;
  unsigned long firstpos, bufsize, windowsize;
  Seqpos currentpos, endpos;
  const Encodedsequence *encseq;
  Encodedsequencescanstate *esr;
};

Windowiterator *windowiterator_new(const Encodedsequence *encseq,
                                   unsigned long windowsize,
                                   Seqpos startpos,
                                   Seqpos endpos)
{
  Windowiterator *wit;

  gt_assert(encseq != NULL);
  gt_assert(endpos <= getencseqtotallength(encseq));
  wit = gt_malloc(sizeof (*wit));
  wit->buffer = gt_malloc(sizeof (Uchar) * windowsize);
  wit->firstpos = wit->bufsize = 0;
  wit->windowsize = windowsize;
  wit->currentpos = startpos;
  wit->endpos = endpos;
  wit->esr = newEncodedsequencescanstate();
  wit->encseq = encseq;
  initEncodedsequencescanstate(wit->esr,encseq,Forwardmode,startpos);
  return wit;
}

void windowiterator_delete(Windowiterator *wit)
{
  gt_free(wit->buffer);
  gt_free(wit);
}

const Uchar *windowiterator_next(Seqpos *currentpos,unsigned long *firstpos,
                                 Windowiterator *wit)
{
  Uchar currentchar;

  while (wit->currentpos < wit->endpos)
  {
    currentchar = sequentialgetencodedchar(wit->encseq,
                                           wit->esr,
                                           wit->currentpos,Forwardmode);
    if (ISSPECIAL(currentchar))
    {
      wit->bufsize = wit->firstpos = 0;
    } else
    {
      if (wit->bufsize < wit->windowsize)
      {
        wit->buffer[wit->bufsize++] = currentchar;
      } else
      {
        wit->buffer[wit->firstpos++] = currentchar;
        if (wit->firstpos == wit->windowsize)
        {
          wit->firstpos = 0;
        }
      }
    }
    if (wit->bufsize == wit->windowsize)
    {
      gt_assert(wit->currentpos >= (Seqpos) (wit->windowsize-1));
      gt_assert(wit->firstpos < wit->windowsize);
      *currentpos = wit->currentpos++;
      *firstpos = wit->firstpos;
      return wit->buffer;
    }
    wit->currentpos++;
  }
  return NULL;
}
