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

#include "core/error.h"
#include "core/fasta_reader.h"
#include "core/seq.h"
#include "core/str.h"

/* GT_Bioseq file endings */
#define GT_BIOSEQ_INDEX        ".gt_bsi"
#define GT_BIOSEQ_RAW          ".gt_bsr"
#define GT_BIOSEQ_FINGERPRINTS ".gt_bsf"

typedef struct GT_Bioseq GT_Bioseq;

/* Construct a new bioseq object (and create the bioseq files, if necessary). */
GT_Bioseq*       gt_bioseq_new(const char *sequence_file, GtError*);
/* Construct a new bioseq object (and always create the the bioseq files). */
GT_Bioseq*       gt_bioseq_new_recreate(const char *sequence_file, GtError*);
GT_Bioseq*       gt_bioseq_new_str(GtStr* sequence_file, GtError*);
/* Construct a new bioseq object (and always create the bioseq files)
   with a certain <fasta_reader>. */
GT_Bioseq*       gt_bioseq_new_with_fasta_reader(const char *sequence_file,
                                           GT_FastaReaderType fasta_reader,
                                           GtError*);
void          gt_bioseq_delete(GT_Bioseq*);
GtAlpha*        gt_bioseq_get_alpha(GT_Bioseq*);
Seq*          gt_bioseq_get_seq(GT_Bioseq*, unsigned long);
const char*   gt_bioseq_get_description(GT_Bioseq*, unsigned long);
/* Return sequence with given <index> (not '\0' terminated). */
const char*   gt_bioseq_get_sequence(GT_Bioseq*, unsigned long index);
const char*   gt_bioseq_get_raw_sequence(GT_Bioseq*);
/* Return MD5 fingerprint of sequence with given <index>. */
const char*   gt_bioseq_get_md5_fingerprint(GT_Bioseq*, unsigned long index);
unsigned long gt_bioseq_get_sequence_length(GT_Bioseq*, unsigned long);
unsigned long gt_bioseq_get_raw_sequence_length(GT_Bioseq*);
unsigned long gt_bioseq_number_of_sequences(GT_Bioseq*);

/* Shows a bioseq on stdout (in fasta format).
   If width is != 0 the sequences are formatted accordingly. */
void gt_bioseq_show_as_fasta(GT_Bioseq*, unsigned long width);

/* Shows a sequence with number ``seqnum'' from a bioseq on stdout (in fasta
   format). If width is != 0 the sequences are formatted accordingly. */
void gt_bioseq_show_sequence_as_fasta(GT_Bioseq*, unsigned long seqnum,
                                   unsigned long width);

/* Shows GC-content on stdout (for DNA files). */
void gt_bioseq_show_gc_content(GT_Bioseq*);

/* Shows bioseq statistics (on stdout). */
void gt_bioseq_show_stat(GT_Bioseq*);

/* Shows bioseq sequence length distribution (on stdout). */
void gt_bioseq_show_seqlengthdistri(GT_Bioseq*);

#endif
