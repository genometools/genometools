/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FASTA_READER_H
#define FASTA_READER_H

#include "str.h"

typedef struct Fasta_reader Fasta_reader;

/* gets called for each description (the start of a fasta entry) */
typedef void (*Fasta_reader_proc_description)(Str*, void *data);
/* gets called for each character of a fasta entry */
typedef void (*Fasta_reader_proc_character)(char, void *data);
/* gets called after a fasta entry has been read */
typedef void (*Fasta_reader_proc_sequence_length)(unsigned long, void *data);

Fasta_reader* fasta_reader_new(Str *sequence_filename);
void          fasta_reader_run(Fasta_reader*, Fasta_reader_proc_description,
                               Fasta_reader_proc_character,
                               Fasta_reader_proc_sequence_length, void *data);
void          fasta_reader_free(Fasta_reader*);

#endif
