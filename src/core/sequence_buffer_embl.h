/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQUENCE_BUFFER_EMBL_H
#define SEQUENCE_BUFFER_EMBL_H

#include "core/sequence_buffer.h"
#include "core/str_array_api.h"

/* implements the ``sequence buffer'' interface for EMBL files */
typedef struct GtSequenceBufferEMBL GtSequenceBufferEMBL;

const GtSequenceBufferClass* gt_sequence_buffer_embl_class(void);
GtSequenceBuffer*            gt_sequence_buffer_embl_new(const GtStrArray*);

bool                         gt_sequence_buffer_embl_guess(const char* txt);

#endif
