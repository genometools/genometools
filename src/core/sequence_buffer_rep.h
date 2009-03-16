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

#ifndef SEQUENCE_BUFFER_REP_H
#define SEQUENCE_BUFFER_REP_H

#include <stdio.h>
#include "core/arraydef.h"
#include "core/genfile.h"
#include "core/queue.h"
#include "core/sequence_buffer.h"
#include "core/str_array.h"

#define INBUFSIZE  4096
#define OUTBUFSIZE 4096

struct GtSequenceBufferClass {
  size_t        size;
  int           (*advance)(GtSequenceBuffer*, GtError*);
  unsigned long (*get_file_index)(GtSequenceBuffer*);
  void          (*free)(GtSequenceBuffer*);
};

typedef struct GtSequenceBufferMembers GtSequenceBufferMembers;

struct GtSequenceBuffer {
  const GtSequenceBufferClass *c_class;
  GtSequenceBufferMembers *pvt;
};

struct GtSequenceBufferMembers {
  unsigned int filenum;
  uint64_t linenum;
  bool complete;
  Filelengthvalues *filelengthtab;
  GtQueue *descptr;
  GtGenFile *inputstream;
  unsigned long reference_count,
                *chardisttab,
                currentfillpos,
                currentinpos,
                nextread,
                nextfree;
  uint64_t lastspeciallength;
  unsigned long long counter;
  const GtStrArray *filenametab;
  Uchar inbuf[INBUFSIZE],
        outbuf[OUTBUFSIZE];
  const Uchar *symbolmap;
};

GtSequenceBuffer* gt_sequence_buffer_create(const GtSequenceBufferClass*);
void*             gt_sequence_buffer_cast(const GtSequenceBufferClass*,
                                          GtSequenceBuffer*);

#endif
