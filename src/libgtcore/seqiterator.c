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

#include "libgtcore/arraydef.h"
#include "libgtcore/chardef.h"
#include "libgtcore/fastabuffer.h"
#include "libgtcore/minmax.h"
#include "libgtcore/strarray.h"
#include "libgtcore/symboldef.h"
#include "libgtcore/seqiterator.h"

struct SeqIterator
{
  FastaBuffer *fb;
  const StrArray *filenametab;
  const Uchar *symbolmap;
  Queue *descptr;
  ArrayUchar sequencebuffer;
  unsigned long long unitnum;
  bool withsequence, exhausted;
  unsigned long long currentread,
                     maxread;
};

SeqIterator* seqiterator_new(const StrArray *filenametab,
                             const Uchar *symbolmap,
                             bool withsequence)
{
  SeqIterator *seqit;
  seqit = ma_malloc(sizeof (SeqIterator));
  INITARRAY(&seqit->sequencebuffer, Uchar);
  seqit->descptr = queue_new();
  seqit->fb = fastabuffer_new(filenametab,
                              symbolmap,
                              false,
                              NULL,
                              seqit->descptr,
                              NULL);
  seqit->exhausted = false;
  seqit->unitnum = 0;
  seqit->withsequence = withsequence;
  seqit->currentread = 0;
  seqit->maxread = 0;
  return seqit;
}

int seqiterator_next(SeqIterator *seqit,
                     const Uchar **sequence,
                     unsigned long *len,
                     char **desc,
                     Error *e)
{
  Uchar charcode;
  int retval;
  bool haserr = false, foundseq = false;

  assert(seqit && len && desc);
  assert((sequence && seqit->withsequence) || !seqit->withsequence);

  if (seqit->exhausted)
  {
    return 0;
  }
  while (true)
  {
    retval = fastabuffer_next(seqit->fb,&charcode,e);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      seqit->exhausted = true;
      break;
    }
    if (seqit->currentread < seqit->maxread)
    {
      seqit->currentread++;
    }
    if (charcode == (Uchar) SEPARATOR)
    {
      if (seqit->sequencebuffer.nextfreeUchar == 0 && seqit->withsequence)
      {
        error_set(e,"sequence %llu is empty", seqit->unitnum);
        haserr = true;
        break;
      }
      *desc = queue_get(seqit->descptr);
      *len = seqit->sequencebuffer.nextfreeUchar;
      if (seqit->withsequence)
      {
        /* make sure the outgoing sequence is '\0' terminated */
        seqit->sequencebuffer.spaceUchar
          [seqit->sequencebuffer.nextfreeUchar] = '\0';
        *sequence = seqit->sequencebuffer.spaceUchar;
      }
      seqit->sequencebuffer.nextfreeUchar = 0;
      foundseq = true;
      seqit->unitnum++;
      break;
    }
    if (seqit->withsequence)
    {
      STOREINARRAY(&seqit->sequencebuffer, Uchar,
                   MAX(1024, seqit->sequencebuffer.nextfreeUchar * 0.5),
                   charcode);
    } else
    {
      seqit->sequencebuffer.nextfreeUchar++;
    }
  }
  if (!haserr && seqit->sequencebuffer.nextfreeUchar > 0)
  {
    *desc = queue_get(seqit->descptr);
    if (seqit->withsequence)
    {
      /* make sure the outgoing sequence is '\0' terminated */
      seqit->sequencebuffer.spaceUchar
        [seqit->sequencebuffer.nextfreeUchar] = '\0';
      *sequence = seqit->sequencebuffer.spaceUchar;
    }
    *len = seqit->sequencebuffer.nextfreeUchar;
    foundseq = true;
    seqit->sequencebuffer.nextfreeUchar = 0;
  }
  if (haserr)
  {
    return -1;
  }
  if (foundseq)
  {
    return 1;
  }
  return 0;
}

const unsigned long long *seqiterator_getcurrentcounter(SeqIterator *seqit,
                                                        unsigned long long
                                                        maxread)
{
  seqit->maxread = maxread;
  return &seqit->currentread;
}

void seqiterator_delete(SeqIterator *seqit)
{
  if (!seqit) return;
  queue_delete_with_contents(seqit->descptr);
  fastabuffer_delete(seqit->fb);
  FREEARRAY(&seqit->sequencebuffer, Uchar);
  seqit->currentread = seqit->maxread;
  ma_free(seqit);
}
