/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef READS2TWOBIT_H
#define READS2TWOBIT_H

/* Reads2Twobit is a fast specialized encoder for large collections of
 * sequencing reads in MultiFasta format. The information can then be
 * accessed as a pointer to a twobitencoding representation in memory
 * or output to disk in a format compatible with GtEncseq.
 *
 * It is designed to be fast and memory efficient, but it is way less flexible
 * than the encoding feature of GtEncseq; in particular by Reads2Twobit:
 * - the input format must be fasta
 * - only nucleic acid sequences are supported
 * - the description lines are discarded
 * - reads containing ambiguities are discarded
 * - there is no support for md5 signatures
 * - supports paired ends libraries: mate pairs must be input as sequences
 *   with the same sequence number in two distinct fasta input files;
 *   they are encoded as subsequent sequences in the encoded sequence
 * - the filetab of the encseq contains "virtual" files (the sequences
 *   are not encoded in the same order as the original files);
 *   in particular each file correspond to a sequencing reads library;
 *   the filename is set to indexname-<libnum>-<insertlength> where <libnum> is
 *   a counter, <insertlength> is the insert length of the library;
 *
 * Both equal lenght and variable length sequence collections are supported.
 * Reads2Twobit automatically switches to variable length mode when the first
 * sequence is encoded, whose length is not the same of the previous sequences.
 */

#include "core/str.h"

typedef struct GtReads2Twobit GtReads2Twobit;

GtReads2Twobit* gt_reads2twobit_new(GtStr *indexname);

void gt_reads2twobit_delete(GtReads2Twobit *r2t);

#define GT_READS2TWOBIT_LIBNAMESEP ':'

/* adds a pair of filenames to the collection, containing reads with pairing
 * information; the two files must contain exactly the same amount of reads;
 * the pair of a read in the first file must be found in the second file at
 * the same sequence number */
void gt_reads2twobit_add_paired(GtReads2Twobit *r2t, const GtStr *filename1,
    const GtStr *filename2, unsigned long insertlength);

/* adds a filename to the collection, containing reads lacking pairing
 * information */
#define gt_reads2twobit_add_unpaired(r2t, filename) \
  gt_reads2twobit_add_paired(r2t, filename, NULL, 0)

/* if libspec is in the form filename1:filename2:insertlength, then parses
 * it and calls add_paired, otherwise calls add_unpaired; returns -1 if
 * the libspec format is not valid */
int gt_reads2twobit_add_library(GtReads2Twobit *r2t, const GtStr *libspec,
                                GtError *err);

/* encodes the sequences in the twobit-encoding format in memory;
 * can be called only once; returns 0 on success */
int gt_reads2twobit_encode(GtReads2Twobit *r2t, GtError *err);

/* (sets the sep characters to the less frequent char if needed, then)
 * writes the sequence collection to disk in GtEncseq format;
 * it must be called after gt_reads2twobit_encode;
 * returns 0 on success */
int gt_reads2twobit_write_encseq(GtReads2Twobit *r2t, GtError *err);

/* writes the sequence collection to disk in MultiFasta format;
 * it must be called after gt_reads2twobit_encode; if <skip> is not NULL,
 * then skips any sequence for which the corresponding bit is set */
int gt_reads2twobit_write_fasta(const GtReads2Twobit *r2t, char *path,
    GtBitsequence *skip, GtError *err);

/* decodes the specified sequence in Fasta format; the <decoded> buffer
 * must be large enough */
void gt_reads2twobit_decode_sequence(const GtReads2Twobit *r2t,
    unsigned long seqnum, char *decoded);

/* decodes the sequences <seqnum_from> to <seqnum_from>+<nofseqs>-1
 * in MultiFasta format and outputs to <outfp>; if <skip> is not NULL,
 * then skips any sequence for which the corresponding bit is set */
void gt_reads2twobit_decode_range(const GtReads2Twobit *r2t,
    GtFile *outfp, unsigned long seqnum_from, unsigned long nofseqs,
    const GtBitsequence *skip);

/* writes the sequence <seqnum> to <outputbuffer>; starts writing at
 * the <outputoffset>-th character encoded by the <outputbuffer> code;
 *
 * returns the number of codes which have been written and sets
 * <*lastcodeoffsetptr> to the offset of the last code
 * (which can be used as <outputoffset> for subsequent calls to the function)
 */
unsigned long gt_reads2twobit_write_encoded(GtReads2Twobit *r2t,
    unsigned long seqnum, GtTwobitencoding *outputbuffer,
    GtTwobitencoding outputoffset, GtTwobitencoding *lastcodeoffsetptr);

/* outputs the seppos array to <path>;
 * if <skip> is not NULL, then skips any position
 * for which the corresponding bit is set */
int gt_reads2twobit_write_seppos(GtReads2Twobit *r2t, char* path,
    GtBitsequence *skip, GtError *err);

/* delete the sequences for which a bit is set in the list;
 * for paired reads: both members of a pair are deleted */
void gt_reads2twobit_delete_sequences(GtReads2Twobit *r2t, GtBitsequence *list);

/* pointer to the internal twobitencoding representation;
 * it must be called after gt_reads2twobit_encode */
GtTwobitencoding *gt_reads2twobit_export_twobitencoding(
    const GtReads2Twobit *r2t);

/* pointer to the internal seppos array;
 * it must be called after gt_reads2twobit_encode */
unsigned long *gt_reads2twobit_export_seppos(const GtReads2Twobit *r2t);

/* encoding statistics, must be called after gt_reads2twobit_encode */
unsigned long gt_reads2twobit_nofseqs(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_seqlen_eqlen(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_seqlen_max(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_seqlen_min(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_total_seqlength(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_nof_invalid_seqs(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_invalid_seqs_totallength(
    const GtReads2Twobit *r2t);

#endif
