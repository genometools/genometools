/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef BIOSEQ_H
#define BIOSEQ_H

#include <libgt/env.h>
#include <libgt/seq.h>
#include <libgt/str.h>

/* Bioseq file endings */
#define GT_BIOSEQ_INDEX ".gt_bsi"
#define GT_BIOSEQ_RAW   ".gt_bsr"

typedef struct Bioseq Bioseq;

/* construct a new bioseq object (and create the bioseq files, if necessary) */
Bioseq*       bioseq_new(const char *sequence_file, Env*);
/* construct a new bioseq object (and always create the the bioseq files) */
Bioseq*       bioseq_new_recreate(const char *sequence_file, Env*);
Bioseq*       bioseq_new_str(Str* sequence_file, Env*);
Alpha*        bioseq_get_alpha(Bioseq*, Env*);
Seq*          bioseq_get_seq(Bioseq*, unsigned long, Env*);
const char*   bioseq_get_description(Bioseq*, unsigned long);
const char*   bioseq_get_sequence(Bioseq*, unsigned long);
const char*   bioseq_get_raw_sequence(Bioseq*);
unsigned long bioseq_get_sequence_length(Bioseq*, unsigned long);
unsigned long bioseq_get_raw_sequence_length(Bioseq*);
unsigned long bioseq_number_of_sequences(Bioseq*);
void          bioseq_delete(Bioseq*, Env*);

/* shows a bioseq on stdout (in fasta format).
   If width is != 0 the sequences are formatted accordingly */
void bioseq_show_as_fasta(Bioseq*, unsigned long width);

/* shows a sequence with number ``seqnum'' from a bioseq on stdout (in fasta
   format). If width is != 0 the sequences are formatted accordingly */
void bioseq_show_sequence_as_fasta(Bioseq*, unsigned long seqnum,
                                   unsigned long width);

/* shows GC-content on stdout (for DNA files */
void bioseq_show_gc_content(Bioseq*, Env*);

/* shows bioseq statistics (on stdout) */
void bioseq_show_stat(Bioseq*);

#endif
