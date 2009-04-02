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

#ifndef SEQITERATOR_H
#define SEQITERATOR_H

#include <inttypes.h>
#include "core/queue.h"
#include "core/str_array.h"
#include "core/sequence_buffer.h"
#include "core/symboldef.h"

typedef struct GtSeqIterator GtSeqIterator;

/* Create a new <GtSeqIterator> for all FASTA files in <filenametab>.
   If a <symbolmap> is given, all read in sequences are transformed with it.
   If <withsequence> equals <true>, FASTA sequences and descriptions are
   processed (otherwise only the descriptions). */
GtSeqIterator* gt_seqiterator_new(const GtStrArray *filenametab,
                                  const Uchar *symbolmap,
                                  bool withsequence);

/* Create a new <GtSeqIterator> for files in <filenametab> using the
   <GtSequenceBuffer> implementation <buffer>.
   If a <symbolmap> is given, all read in sequences are transformed with it.
   If <withsequence> equals <true>, sequences and descriptions are processed
   (otherwise only the descriptions). */
GtSeqIterator* gt_seqiterator_new_with_buffer(GtSequenceBuffer *buffer,
                                              const Uchar *symbolmap,
                                              bool withsequence);

/* Get next <sequence> (of length <len>) and <description> from <seq_iterator>.
   The caller is responsible to free the received <description>.
   Returns 1, if another sequence could be parsed. 0, if all given sequence
   files are exhausted. And -1, if an error occured (<err> is set
   accordingly). */
int            gt_seqiterator_next(GtSeqIterator *seq_iterator,
                                   const Uchar **sequence,
                                   unsigned long *len,
                                   char **description, GtError*);
const unsigned
long long*     gt_seqiterator_getcurrentcounter(GtSeqIterator*,
                                                unsigned long long);
void           gt_seqiterator_delete(GtSeqIterator*);

#endif
