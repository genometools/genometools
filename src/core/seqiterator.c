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

#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/fastabuffer.h"
#include "core/minmax.h"
#include "core/str_array.h"
#include "core/symboldef.h"
#include "core/seqiterator.h"

struct GtSeqIterator
{
  GtFastaBuffer *fb;
  const GtStrArray *filenametab;
  const Uchar *symbolmap;
  GtQueue *descptr;
  ArrayUchar sequencebuffer;
  unsigned long long unitnum;
  bool withsequence, exhausted;
  unsigned long long currentread,
                     maxread;
};

GtSeqIterator* gt_seqiterator_new(const GtStrArray *filenametab,
                                  const Uchar *symbolmap,
                                  bool withsequence)
{
  GtSeqIterator *seqit;
  seqit = gt_malloc(sizeof (GtSeqIterator));
  INITARRAY(&seqit->sequencebuffer, Uchar);
  seqit->descptr = gt_queue_new();
  seqit->fb = gt_fastabuffer_new(filenametab,
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

int gt_seqiterator_next(GtSeqIterator *seqit,
                        const Uchar **sequence,
                        unsigned long *len,
                        char **desc,
                        GtError *err)
{
  Uchar charcode;
  int retval;
  bool haserr = false, foundseq = false;

  gt_assert(seqit && len && desc);
  gt_assert((sequence && seqit->withsequence) || !seqit->withsequence);

  if (seqit->exhausted)
  {
    return 0;
  }
  while (true)
  {
    retval = gt_fastabuffer_next(seqit->fb,&charcode,err);
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
        gt_error_set(err,"sequence %llu is empty", seqit->unitnum);
        haserr = true;
        break;
      }
      *desc = gt_queue_get(seqit->descptr);
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
    *desc = gt_queue_get(seqit->descptr);
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

const unsigned long long *gt_seqiterator_getcurrentcounter(GtSeqIterator *seqit,
                                                           unsigned long long
                                                           maxread)
{
  seqit->maxread = maxread;
  return &seqit->currentread;
}

void gt_seqiterator_delete(GtSeqIterator *seqit)
{
  if (!seqit) return;
  gt_queue_delete_with_contents(seqit->descptr);
  gt_fastabuffer_delete(seqit->fb);
  FREEARRAY(&seqit->sequencebuffer, Uchar);
  seqit->currentread = seqit->maxread;
  gt_free(seqit);
}
