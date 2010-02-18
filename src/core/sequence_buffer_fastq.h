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

#ifndef SEQUENCE_BUFFER_FASTQ_H
#define SEQUENCE_BUFFER_FASTQ_H

#include "core/sequence_buffer.h"
#include "core/str_array_api.h"

/* Implements the ``sequence buffer'' interface for FastQ files.
   TODO: This implementation uses a GtSeqIteratorFastQ, which requires copying
   sequences. It should eventually be replaced by a true FastQ parser/buffer,
   eliminating the need for a separate GtSeqIteratorFastQ (and providing better
   performance). */
typedef struct GtSequenceBufferFastQ GtSequenceBufferFastQ;

const GtSequenceBufferClass* gt_sequence_buffer_fastq_class(void);
GtSequenceBuffer*            gt_sequence_buffer_fastq_new(const GtStrArray*);

bool                         gt_sequence_buffer_fastq_guess(const char* txt);

#endif
