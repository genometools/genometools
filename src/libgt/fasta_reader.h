/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FASTA_READER_H
#define FASTA_READER_H

#include "str.h"

typedef struct FastaReader FastaReader;

/* gets called for each description (the start of a fasta entry) */
typedef void (*FastaReader_proc_description)(Str*, void *data);
/* gets called for each character of a fasta entry */
typedef void (*FastaReader_proc_character)(char, void *data);
/* gets called after a fasta entry has been read */
typedef void (*FastaReader_proc_sequence_length)(unsigned long, void *data);

FastaReader* fasta_reader_new(Str *sequence_filename);
int          fasta_reader_run(FastaReader*, FastaReader_proc_description,
                              FastaReader_proc_character,
                              FastaReader_proc_sequence_length, void *data,
                              Env*);
void         fasta_reader_delete(FastaReader*, Env*);

#endif
