/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <libgtcore/str.h>

typedef struct FastaReader FastaReader;

/* gets called for each description (the start of a fasta entry) */
typedef int (*FastaReaderProcDescription)(Str*, void *data, Env*);
/* gets called for each character of a fasta entry */
typedef int (*FastaReaderProcCharacter)(char, void *data, Env*);
/* gets called after a fasta entry has been read */
typedef int (*FastaReaderProcSequenceLength)(unsigned long, void *data, Env*);

/* construct a new fasta reader for the file named <sequence_filename>, pass
   NULL to read from stdin */
FastaReader* fasta_reader_new(Str *sequence_filename, Env*);
int          fasta_reader_run(FastaReader*, FastaReaderProcDescription,
                              FastaReaderProcCharacter,
                              FastaReaderProcSequenceLength, void *data, Env*);
void         fasta_reader_delete(FastaReader*, Env*);

#endif
