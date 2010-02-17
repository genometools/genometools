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

#ifndef SEQITERATOR_SEQUENCE_BUFFER_H
#define SEQITERATOR_SEQUENCE_BUFFER_H

#include "core/error_api.h"
#include "core/seqiterator.h"
#include "core/str_array.h"
#include "core/sequence_buffer.h"

typedef struct GtSeqIteratorSequenceBuffer GtSeqIteratorSequenceBuffer;

/* Create a new <GtSeqIterator> for all sequence files in <filenametab>.
   All files have to be of the same format, which will be guessed by examining
   the beginning of the first file. If an error occurs, NULL is returned (see
   the <err> object for details). */
GtSeqIterator* gt_seqiterator_sequence_buffer_new(const GtStrArray *filenametab,
                                                  GtError *err);

/* Create a new <GtSeqIterator> for files in <filenametab> using the
   <GtSequenceBuffer> implementation <buffer>. */
GtSeqIterator* gt_seqiterator_sequence_buffer_new_with_buffer(
                                                      GtSequenceBuffer *buffer);

#endif
