/*
   Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
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

#ifndef BIOSEQ_API_H
#define BIOSEQ_API_H

#include "core/alphabet_api.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/fasta_reader_api.h"
#include "core/seq_api.h"
#include "core/str_api.h"

/* <GtBioseq> represents a simple collection of biosequences. */
typedef struct GtBioseq GtBioseq;

/* Construct a new <GtBioseq> object (and create the bioseq files, if
   necessary). */
GtBioseq*   gt_bioseq_new(const char *sequence_file, GtError*);
/* Construct a new <GtBioseq> object (and always create the the bioseq
   files). */
GtBioseq*   gt_bioseq_new_recreate(const char *sequence_file, GtError*);
/* Construct a new <GtBioseq> object (and always create the the bioseq
   files, if necessary). Filename is given as a <GtString>. */
GtBioseq*   gt_bioseq_new_str(GtStr* sequence_file, GtError*);
/* Delete the <bioseq>. */
void        gt_bioseq_delete(GtBioseq *bioseq);
/* Delete the index files belonging to <bioseq>. */
void        gt_bioseq_delete_indices(GtBioseq *bioseq);
/* Return the <GtAlphabet> associated with <bioseq>. */
GtAlphabet* gt_bioseq_get_alphabet(GtBioseq *bioseq);
/* Return sequence with given <index> as <GtSeq>. */
GtSeq*      gt_bioseq_get_seq(GtBioseq *bioseq, GtUword index);
/* Return subsequence within the boundaries <start> and <end> of the sequence
   with given <index> as <GtSeq>. */
GtSeq*      gt_bioseq_get_seq_range(GtBioseq*, GtUword index,
                                    GtUword start, GtUword end);
/* Return description of the sequence with given <index>. */
const char* gt_bioseq_get_description(GtBioseq*, GtUword);
/* Return character at position <position> of the sequence with given
   <index>. */
char        gt_bioseq_get_char(const GtBioseq*, GtUword index,
                               GtUword position);
/* Return TRUE if sequence with given <index> contains wildcards. */
bool        gt_bioseq_seq_has_wildcards(const GtBioseq *bioseq,
                                        GtUword idx);
/* Return sequence with given <index> (not '\0' terminated). */
char*       gt_bioseq_get_sequence(const GtBioseq *bioseq, GtUword index);
/* Return subsequence within the boundaries <start> and <end> of the sequence
   with given <index> (not '\0' terminated). */
char*       gt_bioseq_get_sequence_range(const GtBioseq *bioseq, GtUword index,
                                         GtUword start,
                                         GtUword end);
/* Return character at position <position> of the sequence with given
   <index>. Character is encoded according to alphabet. */
GtUchar     gt_bioseq_get_encoded_char(const GtBioseq *bioseq, GtUword index,
                                       GtUword position);
/* Writes encoded sequence with given <index> as to <out>. The sequence
   is encoded according to alphabet. */
void        gt_bioseq_get_encoded_sequence(const GtBioseq *bioseq,
                                           GtUchar *out,
                                           GtUword index);
/* Return subsequence within the boundaries <start> and <end> of the sequence
   with given <index> (not '\0' terminated). The sequence is encoded according
   to the alphabet of <bioseq>. */
void        gt_bioseq_get_encoded_sequence_range(const GtBioseq *bioseq,
                                                 GtUchar *out,
                                                 GtUword index,
                                                 GtUword start,
                                                 GtUword end);
/* Return MD5 fingerprint of sequence with given <index>. */
const char* gt_bioseq_get_md5_fingerprint(GtBioseq *bioseq, GtUword index);
/* Return filename of sequence file underlying <bioseq>. */
const char* gt_bioseq_filename(const GtBioseq *bioseq);
/* Return length of the sequence with given <index> in <bioseq>. */
GtUword     gt_bioseq_get_sequence_length(const GtBioseq *bioseq,
                                          GtUword index);
/* Return length of all sequences in <bioseq>. */
GtUword     gt_bioseq_get_total_length(const GtBioseq *bioseq);
/* Return count of all sequences in <bioseq>. */
GtUword     gt_bioseq_number_of_sequences(GtBioseq *bioseq);
/* Return the index of the (first) sequence with given <MD5> contained in
   <bioseq>, if it exists. Otherwise <GT_UNDEF_UWORD> is returned. */
GtUword     gt_bioseq_md5_to_index(GtBioseq *bioseq, const char *MD5);
/* Shows a <bioseq> on <outfp> (in fasta format).
   If <width> is != 0 the sequences are formatted accordingly. */
void        gt_bioseq_show_as_fasta(GtBioseq *bioseq, GtUword width,
                                    GtFile *outfp);
/* Shows a sequence with number <seqnum> from a <bioseq> on <outfp> (in fasta
   format). If <width> is != 0 the sequences are formatted accordingly. */
void        gt_bioseq_show_sequence_as_fasta(GtBioseq *bioseq, GtUword seqnum,
                                             GtUword width, GtFile *outfp);
/* Shows GC-content on <outfp> (for DNA files). */
void        gt_bioseq_show_gc_content(GtBioseq*, GtFile *outfp);
/* Shows <bioseq> statistics (on <outfp>). */
void        gt_bioseq_show_stat(GtBioseq *bioseq, GtFile *outfp);
/* Shows <bioseq> sequence length distribution (on <outfp>). */
void        gt_bioseq_show_seqlengthdistri(GtBioseq *bioseq, GtFile *outfp);

#endif
