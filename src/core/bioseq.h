/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/alphabet.h"
#include "core/encseq.h"
#include "core/error.h"
#include "core/fasta_reader.h"
#include "core/seq.h"
#include "core/str.h"

typedef struct GtBioseq GtBioseq;

/* Construct a new bioseq object (and create the bioseq files, if necessary). */
GtBioseq*     gt_bioseq_new(const char *sequence_file, GtError*);
/* Construct a new bioseq object (and always create the the bioseq files). */
GtBioseq*     gt_bioseq_new_recreate(const char *sequence_file, GtError*);
GtBioseq*     gt_bioseq_new_str(GtStr* sequence_file, GtError*);
void          gt_bioseq_delete(GtBioseq*);
void          gt_bioseq_delete_indices(GtBioseq*);
GtAlphabet*   gt_bioseq_get_alphabet(GtBioseq*);
GtSeq*        gt_bioseq_get_seq(GtBioseq*, unsigned long);
GtSeq*        gt_bioseq_get_seq_range(GtBioseq*, unsigned long index,
                                      unsigned long start, unsigned long end);
const char*   gt_bioseq_get_description(GtBioseq*, unsigned long);
char          gt_bioseq_get_char(const GtBioseq*, unsigned long index,
                                 unsigned long position);
/* Return sequence with given <index> (not '\0' terminated). */
char*         gt_bioseq_get_sequence(const GtBioseq*, unsigned long index);
char*         gt_bioseq_get_sequence_range(const GtBioseq*, unsigned long index,
                                           unsigned long start,
                                           unsigned long end);
GtUchar       gt_bioseq_get_encoded_char(const GtBioseq*, unsigned long index,
                                         unsigned long position);
void          gt_bioseq_get_encoded_sequence(const GtBioseq*, GtUchar *out,
                                             unsigned long index);
void          gt_bioseq_get_encoded_sequence_range(const GtBioseq*,
                                                   GtUchar *out,
                                                   unsigned long index,
                                                   unsigned long start,
                                                   unsigned long end);
/* Return MD5 fingerprint of sequence with given <index>. */
const char*   gt_bioseq_get_md5_fingerprint(GtBioseq*, unsigned long index);
/* Return filename of sequence file underlying <bioseq>. */
const char*   gt_bioseq_filename(const GtBioseq *bioseq);
unsigned long gt_bioseq_get_sequence_length(const GtBioseq*,
                                            unsigned long index);
unsigned long gt_bioseq_get_total_length(const GtBioseq*);
unsigned long gt_bioseq_number_of_sequences(GtBioseq*);
/* Return the index of the (first) sequence with given <MD5> contained in
   <bioseq>, if it exists. Otherwise <GT_UNDEF_ULONG> is returned. */
unsigned long gt_bioseq_md5_to_index(GtBioseq *bioseq, const char *MD5);

/* Shows a <bioseq> on <outfp> (in fasta format).
   If <width> is != 0 the sequences are formatted accordingly. */
void gt_bioseq_show_as_fasta(GtBioseq *bioseq, unsigned long width,
                             GtFile *outfp);

/* Shows a sequence with number <seqnum> from a <bioseq> on <outfp> (in fasta
   format). If <width> is != 0 the sequences are formatted accordingly. */
void gt_bioseq_show_sequence_as_fasta(GtBioseq *bioseq, unsigned long seqnum,
                                      unsigned long width, GtFile *outfp);

/* Shows GC-content on <outfp> (for DNA files). */
void gt_bioseq_show_gc_content(GtBioseq*, GtFile *outfp);

/* Shows <bioseq> statistics (on <outfp>). */
void gt_bioseq_show_stat(GtBioseq *bioseq, GtFile *outfp);

/* Shows <bioseq> sequence length distribution (on <outfp>). */
void gt_bioseq_show_seqlengthdistri(GtBioseq *bioseq, GtFile *outfp);

#endif
