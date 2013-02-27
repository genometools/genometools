/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/cstr_api.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/sequence_buffer_inline.h"
#include "core/sequence_buffer_rep.h"
#include "core/sequence_buffer_fastq.h"
#include "core/seq_iterator_fastq_api.h"
#include "core/undef_api.h"

#define FASTQ_START_SYMBOL '@'

struct GtSequenceBufferFastQ
{
  const GtSequenceBuffer parent_instance;
  GtSeqIteratorFastQ *seqit;
  const GtStrArray *sequences;
  bool carryseparator;
  GtStr *overflowbuffer;
};

#define gt_sequence_buffer_fastq_cast(SB)\
        gt_sequence_buffer_cast(gt_sequence_buffer_fastq_class(), SB)

static int gt_sequence_buffer_fastq_advance(GtSequenceBuffer *sb, GtError *err)
{
  unsigned long currentoutpos = 0, currentfileadd = 0, currentfileread = 0,
                seqlen = 0, desclen, cnt, newfilenum;
  GtSequenceBufferMembers *pvt;
  GtSequenceBufferFastQ *sbfq;
  const GtUchar *seq;
  char *desc;

  int had_err = 0;

  gt_error_check(err);

  sbfq = gt_sequence_buffer_fastq_cast(sb);
  pvt = sb->pvt;

  if (!sbfq->seqit) {
    sbfq->seqit = (GtSeqIteratorFastQ*)
                         gt_seq_iterator_fastq_new(sbfq->sequences, err);
    if (!sbfq->seqit) return -1;
  }

  /* did the last buffer end with a sequence boundary?
     if so, we need to provide an additional separator! */
  if (sbfq->carryseparator) {
    pvt->outbuf[currentoutpos++] = (GtUchar) SEPARATOR;
    currentfileread++;
    pvt->lastspeciallength++;
    currentfileadd++;
    sbfq->carryseparator = false;
  }

  if (gt_str_length(sbfq->overflowbuffer) > 0) {
    /* we still have surplus sequence from the last line, process that first */
    const char *overflowedstring;
    char cc;
    overflowedstring = gt_str_get(sbfq->overflowbuffer);
    while (!had_err && (cc = *(overflowedstring++)) != '\0')
    {
      if ((had_err = process_char(sb, currentoutpos, cc, err)))
          return had_err;
      currentoutpos++;
      currentfileadd++;
      currentfileread++;
    }
    pvt->outbuf[currentoutpos++] = (GtUchar) SEPARATOR;
    currentfileread++;
    pvt->lastspeciallength++;
    gt_str_reset(sbfq->overflowbuffer);
    gt_assert(gt_str_length(sbfq->overflowbuffer) == 0);
  }

  /* iterate over sequences */
  while (true) {
    had_err = gt_seq_iterator_next((GtSeqIterator*) sbfq->seqit, &seq, &seqlen,
                                  &desc, err);
    newfilenum = gt_seq_iterator_fastq_get_file_index(sbfq->seqit);

    /* update current file info */
    if (pvt->filenum != newfilenum) {
      if (pvt->filelengthtab) {
        pvt->filelengthtab[pvt->filenum].length += (uint64_t) currentfileread;
        pvt->filelengthtab[pvt->filenum].effectivelength += (uint64_t)
                                                               currentfileadd;
      }
      currentfileread = currentfileadd = 0;
      pvt->filenum = newfilenum;
    }

    /* error check! */
    if (had_err < 0) {
      return had_err;
    }

    /* end of iterator reached */
    if (had_err == 0) {
      pvt->complete = true;
      /* remove last separator */
      pvt->outbuf[--currentoutpos] = (GtUchar) '\0';
      currentfileadd--;
      break;
    }

    /* copy sequence */
    for (cnt=0;cnt<seqlen;cnt++) {
      if (currentoutpos >= (unsigned long) OUTBUFSIZE) {
        gt_str_append_char(sbfq->overflowbuffer, seq[cnt]);
      } else {
        if ((had_err = process_char(sb, currentoutpos, seq[cnt], err)))
          return had_err;
        currentoutpos++;
        currentfileadd++;
        currentfileread++;
      }
    }

    /* place separator after sequence (or defer) */
    if (gt_str_length(sbfq->overflowbuffer) == 0) {
      if (currentoutpos >= (unsigned long) OUTBUFSIZE)
        sbfq->carryseparator = true;
      else {
        pvt->outbuf[currentoutpos++] = (GtUchar) SEPARATOR;
        pvt->lastspeciallength++;
        currentfileadd++;
      }
    }

    desclen = strlen(desc)+1;
    currentfileread += desclen; /* XXX: we cannot know whether the desc was */
    pvt->counter += desclen;    /*      repeated before the qualities line  */

    /* enqueue description */
    if (pvt->descptr) {
      int i;
      for (i = 0; i < strlen(desc); i++) {
        gt_desc_buffer_append_char(pvt->descptr, desc[i]);
      }
      gt_desc_buffer_finish(pvt->descptr);
    }

    /* if buffer is full, return it */
    if (currentoutpos >= (unsigned long) OUTBUFSIZE) {
      pvt->nextfree = MIN(currentoutpos, OUTBUFSIZE);
      had_err = 0;
      break;
    }
  }

  /* update information about read characters */
  if (pvt->filelengthtab) {
    pvt->filelengthtab[pvt->filenum].length
      += (uint64_t) currentfileread;
    pvt->filelengthtab[pvt->filenum].effectivelength
      += (uint64_t) currentfileadd;
  }
  pvt->nextfree = MIN(currentoutpos, OUTBUFSIZE);
  return had_err;
}

static void gt_sequence_buffer_fastq_free(GtSequenceBuffer *sb)
{
  GtSequenceBufferFastQ *sbfq = gt_sequence_buffer_fastq_cast(sb);
  gt_seq_iterator_delete((GtSeqIterator*) sbfq->seqit);
  gt_str_delete(sbfq->overflowbuffer);
}

static unsigned long
gt_sequence_buffer_fastq_get_file_index(GtSequenceBuffer *sb)
{
  GtSequenceBufferFastQ *sbfq;
  gt_assert(sb);
  sbfq = gt_sequence_buffer_fastq_cast(sb);
  return gt_seq_iterator_fastq_get_file_index((GtSeqIteratorFastQ*)
                                              sbfq->seqit);
}

const GtSequenceBufferClass* gt_sequence_buffer_fastq_class(void)
{
  static const GtSequenceBufferClass sbc = { sizeof (GtSequenceBufferFastQ),
                                        gt_sequence_buffer_fastq_advance,
                                        gt_sequence_buffer_fastq_get_file_index,
                                        gt_sequence_buffer_fastq_free };
  return &sbc;
}

bool gt_sequence_buffer_fastq_guess(const char* txt)
{
  return (txt[0] == FASTQ_START_SYMBOL);
}

GtSequenceBuffer* gt_sequence_buffer_fastq_new(const GtStrArray *sequences)
{
  GtSequenceBuffer *sb;
  GtSequenceBufferFastQ *sbfq;
  sb = gt_sequence_buffer_create(gt_sequence_buffer_fastq_class());
  sbfq = gt_sequence_buffer_fastq_cast(sb);
  sbfq->seqit = NULL;
  sbfq->sequences = sequences;
  sb->pvt->filenametab = sequences;
  sbfq->overflowbuffer = gt_str_new();
  sb->pvt->filenum = 0;
  sb->pvt->nextread = sb->pvt->nextfree = 0;
  sb->pvt->lastspeciallength = 0;
  return sb;
}
