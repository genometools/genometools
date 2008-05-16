/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef BIOSEQ_H
#define BIOSEQ_H

#include "libgtcore/error.h"
#include "libgtcore/fasta_reader.h"
#include "libgtcore/seq.h"
#include "libgtcore/str.h"

/* Bioseq file endings */
#define GT_BIOSEQ_INDEX        ".gt_bsi"
#define GT_BIOSEQ_RAW          ".gt_bsr"
#define GT_BIOSEQ_FINGERPRINTS ".gt_bsf"

typedef struct Bioseq Bioseq;

/* Construct a new bioseq object (and create the bioseq files, if necessary). */
Bioseq*       bioseq_new(const char *sequence_file, Error*);
/* Construct a new bioseq object (and always create the the bioseq files). */
Bioseq*       bioseq_new_recreate(const char *sequence_file, Error*);
Bioseq*       bioseq_new_str(Str* sequence_file, Error*);
/* Construct a new bioseq object (and always create the bioseq files)
   with a certain <fasta_reader>. */
Bioseq*       bioseq_new_with_fasta_reader(const char *sequence_file,
                                           FastaReaderType fasta_reader,
                                           Error*);
void          bioseq_delete(Bioseq*);
Alpha*        bioseq_get_alpha(Bioseq*);
Seq*          bioseq_get_seq(Bioseq*, unsigned long);
const char*   bioseq_get_description(Bioseq*, unsigned long);
/* Return sequence with given <index> (not '\0' terminated). */
const char*   bioseq_get_sequence(Bioseq*, unsigned long index);
const char*   bioseq_get_raw_sequence(Bioseq*);
/* Return MD5 fingerprint of sequence with given <index>. */
const char*   bioseq_get_md5_fingerprint(Bioseq*, unsigned long index);
unsigned long bioseq_get_sequence_length(Bioseq*, unsigned long);
unsigned long bioseq_get_raw_sequence_length(Bioseq*);
unsigned long bioseq_number_of_sequences(Bioseq*);

/* Shows a bioseq on stdout (in fasta format).
   If width is != 0 the sequences are formatted accordingly. */
void bioseq_show_as_fasta(Bioseq*, unsigned long width);

/* Shows a sequence with number ``seqnum'' from a bioseq on stdout (in fasta
   format). If width is != 0 the sequences are formatted accordingly. */
void bioseq_show_sequence_as_fasta(Bioseq*, unsigned long seqnum,
                                   unsigned long width);

/* Shows GC-content on stdout (for DNA files). */
void bioseq_show_gc_content(Bioseq*);

/* Shows bioseq statistics (on stdout). */
void bioseq_show_stat(Bioseq*);

/* Shows bioseq sequence length distribution (on stdout). */
void bioseq_show_seqlengthdistri(Bioseq*);

#endif
