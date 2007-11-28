/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

/* this class is deprecated, use the SeqIterator class instead! */

#include "libgtcore/str.h"

typedef struct FastaReader FastaReader;

/* gets called for each description (the start of a fasta entry) */
typedef int (*FastaReaderProcDescription)(Str*, void *data, Error*);
/* gets called for each character of a fasta entry */
typedef int (*FastaReaderProcCharacter)(char, void *data, Error*);
/* gets called after a fasta entry has been read */
typedef int (*FastaReaderProcSequenceLength)(unsigned long, void *data, Error*);

/* construct a new fasta reader for the file named <sequence_filename>, pass
   NULL to read from stdin */
FastaReader* fasta_reader_new(Str *sequence_filename);
int          fasta_reader_run(FastaReader*, FastaReaderProcDescription,
                              FastaReaderProcCharacter,
                              FastaReaderProcSequenceLength, void *data,
                              Error*);
void         fasta_reader_delete(FastaReader*);

#endif
