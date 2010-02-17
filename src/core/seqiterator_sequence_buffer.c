/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/minmax.h"
#include "core/seqiterator_rep.h"
#include "core/seqiterator_sequence_buffer.h"

struct GtSeqIteratorSequenceBuffer
{
  const GtSeqIterator parent_instance;
  GtSequenceBuffer *fb;
  const GtStrArray *filenametab;
  const GtUchar *symbolmap;
  GtQueue *descptr;
  GtArrayGtUchar sequencebuffer;
  unsigned long long unitnum;
  bool withsequence, exhausted;
  unsigned long long currentread,
                     maxread;
};

#define gt_seqiterator_sequence_buffer_cast(SI)\
        gt_seqiterator_cast(gt_seqiterator_sequence_buffer_class(), SI);

const GtSeqIteratorClass* gt_seqiterator_sequence_buffer_class(void);

GtSeqIterator* gt_seqiterator_sequence_buffer_new(const GtStrArray *filenametab,
                                                  GtError *err)
{
  GtSeqIterator *si;
  GtSequenceBuffer *sb = gt_sequence_buffer_new_guess_type(filenametab, err);
  if (!sb)
    return NULL;
  si = gt_seqiterator_sequence_buffer_new_with_buffer(sb);
  gt_sequence_buffer_delete(sb); /* drop this reference */
  return si;
}

GtSeqIterator* gt_seqiterator_sequence_buffer_new_with_buffer(
                                                          GtSequenceBuffer *buf)
{
  GtSeqIterator *si;
  GtSeqIteratorSequenceBuffer *seqit;
  si = gt_seqiterator_create(gt_seqiterator_sequence_buffer_class());
  seqit = gt_seqiterator_sequence_buffer_cast(si);
  GT_INITARRAY(&seqit->sequencebuffer, GtUchar);
  seqit->descptr = gt_queue_new();
  seqit->fb = gt_sequence_buffer_ref(buf);
  gt_sequence_buffer_set_desc_queue(seqit->fb, seqit->descptr);
  seqit->exhausted = false;
  seqit->unitnum = 0;
  seqit->withsequence = true;
  seqit->currentread = 0;
  seqit->maxread = 0;
  return si;
}

static void gt_seqiterator_sequence_buffer_set_symbolmap(GtSeqIterator *si,
                                                         const GtUchar
                                                           *symbolmap)
{
  GtSeqIteratorSequenceBuffer *seqit;
  gt_assert(si);
  seqit = gt_seqiterator_sequence_buffer_cast(si);
  gt_sequence_buffer_set_symbolmap(seqit->fb, symbolmap);
}

static void gt_seqiterator_sequence_buffer_set_sequence_output(
                                                          GtSeqIterator *si,
                                                          bool withsequence)
{
  GtSeqIteratorSequenceBuffer *seqit;
  gt_assert(si);
  seqit = gt_seqiterator_sequence_buffer_cast(si);
  seqit->withsequence = withsequence;
}

static int gt_seqiterator_sequence_buffer_next(GtSeqIterator *si,
                                               const GtUchar **sequence,
                                               unsigned long *len,
                                               char **desc,
                                               GtError *err)
{
  GtSeqIteratorSequenceBuffer *seqit;
  GtUchar charcode;
  int retval;
  bool haserr = false, foundseq = false;
  gt_assert(si);
  gt_assert(len && desc);

  seqit = gt_seqiterator_sequence_buffer_cast(si);
  gt_assert((sequence && seqit->withsequence) || !seqit->withsequence);

  if (seqit->exhausted)
  {
    return 0;
  }
  while (true)
  {
    retval = gt_sequence_buffer_next(seqit->fb,&charcode,err);
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
    if (charcode == (GtUchar) SEPARATOR)
    {
      if (seqit->sequencebuffer.nextfreeGtUchar == 0 && seqit->withsequence)
      {
        gt_error_set(err,"sequence %llu is empty", seqit->unitnum);
        haserr = true;
        break;
      }
      *desc = gt_queue_get(seqit->descptr);
      *len = seqit->sequencebuffer.nextfreeGtUchar;
      if (seqit->withsequence)
      {
        /* make sure the outgoing sequence is '\0' terminated */
        seqit->sequencebuffer.spaceGtUchar
          [seqit->sequencebuffer.nextfreeGtUchar] = '\0';
        *sequence = seqit->sequencebuffer.spaceGtUchar;
      }
      seqit->sequencebuffer.nextfreeGtUchar = 0;
      foundseq = true;
      seqit->unitnum++;
      break;
    }
    if (seqit->withsequence)
    {
      GT_STOREINARRAY(&seqit->sequencebuffer, GtUchar,
                   MAX(1024, seqit->sequencebuffer.nextfreeGtUchar * 0.5),
                   charcode);
    } else
    {
      seqit->sequencebuffer.nextfreeGtUchar++;
    }
  }
  if (!haserr && seqit->sequencebuffer.nextfreeGtUchar > 0)
  {
    *desc = gt_queue_get(seqit->descptr);
    if (seqit->withsequence)
    {
      /* make sure the outgoing sequence is '\0' terminated */
      seqit->sequencebuffer.spaceGtUchar
        [seqit->sequencebuffer.nextfreeGtUchar] = '\0';
      *sequence = seqit->sequencebuffer.spaceGtUchar;
    }
    *len = seqit->sequencebuffer.nextfreeGtUchar;
    foundseq = true;
    seqit->sequencebuffer.nextfreeGtUchar = 0;
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

static const unsigned long long*
gt_seqiterator_sequence_buffer_getcurrentcounter(GtSeqIterator *si,
                                                 unsigned long long
                                                 maxread)
{
  GtSeqIteratorSequenceBuffer *seqit;
  gt_assert(si);
  seqit = gt_seqiterator_sequence_buffer_cast(si);
  seqit->maxread = maxread;
  return &seqit->currentread;
}

static void gt_seqiterator_sequence_buffer_delete(GtSeqIterator *si)
{
  if (!si) return;
  GtSeqIteratorSequenceBuffer *seqit;
  gt_assert(si);
  seqit = gt_seqiterator_sequence_buffer_cast(si);
  gt_queue_delete_with_contents(seqit->descptr);
  gt_sequence_buffer_delete(seqit->fb);
  GT_FREEARRAY(&seqit->sequencebuffer, GtUchar);
  seqit->currentread = seqit->maxread;
}

const GtSeqIteratorClass* gt_seqiterator_sequence_buffer_class(void)
{
  static const GtSeqIteratorClass *sic = NULL;
  if (!sic) {
    sic = gt_seqiterator_class_new(sizeof (GtSeqIteratorSequenceBuffer),
                             gt_seqiterator_sequence_buffer_set_symbolmap,
                             gt_seqiterator_sequence_buffer_set_sequence_output,
                             gt_seqiterator_sequence_buffer_next,
                             gt_seqiterator_sequence_buffer_getcurrentcounter,
                             NULL,
                             gt_seqiterator_sequence_buffer_delete);
  }
  return sic;
}
