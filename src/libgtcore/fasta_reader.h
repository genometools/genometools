/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <stdbool.h>
#include "libgtcore/error.h"

/* the ``fasta reader'' interface */
typedef struct FastaReaderClass FastaReaderClass;
typedef struct FastaReader FastaReader;

typedef enum {
  FASTA_READER_REC,
  FASTA_READER_FSM,
  FASTA_READER_SEQIT
} FastaReaderType;

/* gets called for each description (the start of a fasta entry) */
typedef int (*FastaReaderProcDescription)(const char *description,
                                          unsigned long length, void *data,
                                          Error*);
/* gets called for each sequence part of a fasta entry */
typedef int (*FastaReaderProcSequencePart)(const char *seqpart,
                                           unsigned long length, void *data,
                                           Error*);
/* gets called after a fasta entry has been read */
typedef int (*FastaReaderProcSequenceLength)(unsigned long, void *data, Error*);

/* construct a new fasta reader for the file named <sequence_filename>, pass
   NULL to read from stdin */
int          fasta_reader_run(FastaReader*, FastaReaderProcDescription,
                              FastaReaderProcSequencePart,
                              FastaReaderProcSequenceLength, void *data,
                              Error*);
void         fasta_reader_delete(FastaReader*);

#endif
