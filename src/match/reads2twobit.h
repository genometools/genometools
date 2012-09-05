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

#include "core/fptr_api.h"
#include "core/intbits.h"
#include "core/str.h"

/* The <GtReads2Twobit> class is a specialized encoder for large collections of
   sequencing reads in FastQ and Fasta format. After encoding, the information
   is accessible as a GtTwobitencoding representation in memory and can be
   output to disk in a GtEncseq-compatible format.
   It is designed to be fast and memory efficient, and thus has some
   limitations compared to the encoding feature of GtEncseq: the only supported
   input is Fasta/FastQ with a DNA alphabet; reads containing wildchars are
   removed; the descriptions are discarded; the md5 are not calculated.
   Paired-end libraries are supported. The input may consists in sets of two
   files (forward reads and reverse reads, in the same order) or a single file
   with interleaved reads. In paired-end mode, if a sequence is discarded (e.g.
   because it contains ambiguities), then its mate is discarded too.
   Both equal lenght and variable length sequence collections are supported.
   Reads2Twobit automatically switches to variable length mode when the first
   sequence is encoded, whose length is not the same of the previous
   sequences. */

typedef struct GtReads2Twobit GtReads2Twobit;

GtReads2Twobit* gt_reads2twobit_new(GtStr *indexname);

void gt_reads2twobit_delete(GtReads2Twobit *r2t);

#define GT_READS2TWOBIT_LIBSPECSEP ':'
#define GT_READS2TWOBIT_INSERTSEP  ','

/* Adds the library described by <libspec> to the collection.
   If <libspec> contains at least a ':' is is assumed to be in the
   form "filename1[:filename2]:insertlength[,stdev]" and parsed;
   on success 0 is returned, if a parsing error occurs, -1 is returned.
   If <libspec> does not contain any ':', it is assumed to be a filename
   of a single-end library; in this case the function returns 0. */
int gt_reads2twobit_add_library(GtReads2Twobit *r2t,
                                const GtStr *libspec,
                                GtError *err);

/* Use phred64 scores instead of phred33;
   it must be called before <gt_reads2twobit_encode>. */
void gt_reads2twobit_use_phred64(GtReads2Twobit *r2t);

/* filter those reads which contain more than <maxlow> positions
   whose quality is no more than <lowqual>. */
void gt_reads2twobit_set_quality_filter(GtReads2Twobit *r2t,
    unsigned long maxlow, char lowqual);

/* Encodes the sequences in the twobit-encoding format in memory;
   can be called only once; returns 0 on success, a negative number on
   error and sets <err> accordingly. */
int gt_reads2twobit_encode(GtReads2Twobit *r2t, GtError *err);

/* encoding statistics, must be called after <gt_reads2twobit_encode> */
unsigned long gt_reads2twobit_nofseqs(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_seqlen_eqlen(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_seqlen_max(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_seqlen_min(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_total_seqlength(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_nof_invalid_seqs(const GtReads2Twobit *r2t);
unsigned long gt_reads2twobit_invalid_seqs_totallength(
    const GtReads2Twobit *r2t);

/* Writes the sequence collection to disk in a GtEncseq-compatible format;
   sets the separator positions to the less frequent character when needed;
   it must be called after <gt_reads2twobit_encode>; returns 0 on success. */
int gt_reads2twobit_write_encseq(GtReads2Twobit *r2t, GtError *err);

/* writes the sequence collection to disk in MultiFasta format;
   it must be called after <gt_reads2twobit_encode>; if <skip> is not NULL,
   then skips any sequence for which the corresponding bit is set */
int gt_reads2twobit_write_fasta(const GtReads2Twobit *r2t, char *path,
    GtBitsequence *skip, GtError *err);

/* decodes the specified sequence in Fasta format; the <decoded> buffer
   must be large enough */
void gt_reads2twobit_decode_sequence(const GtReads2Twobit *r2t,
    unsigned long seqnum, char *decoded);

/* decodes the sequences <seqnum_from> to <seqnum_from>+<nofseqs>-1
   in MultiFasta format and outputs to <outfp>; if <skip> is not NULL,
   then skips any sequence for which the corresponding bit is set */
void gt_reads2twobit_decode_range(const GtReads2Twobit *r2t,
    GtFile *outfp, unsigned long seqnum_from, unsigned long nofseqs,
    const GtBitsequence *skip);

/* writes the sequence <seqnum> to <outputbuffer>; starts writing at
  the <outputoffset>-th character encoded by the <outputbuffer> code;
  returns a pointer to the next buffer position where a code is not complete
  and sets <*lastcodeoffsetptr> to the offset of the last code
  (which can be used as <outputoffset> for subsequent calls to the function)
 */
GtTwobitencoding* gt_reads2twobit_write_encoded(GtReads2Twobit *r2t,
    unsigned long seqnum, GtTwobitencoding *outputbuffer,
    GtTwobitencoding outputoffset, GtTwobitencoding *lastcodeoffsetptr);

/* outputs the seppos array to <path>;
  if <skip> is not NULL, then skips any position
  for which the corresponding bit is set */
int gt_reads2twobit_write_seppos(GtReads2Twobit *r2t, char* path,
    GtBitsequence *skip, GtError *err);

/* delete the sequences for which a bit is set in the list;
   for paired reads: both members of a pair are deleted */
void gt_reads2twobit_delete_sequences(GtReads2Twobit *r2t, GtBitsequence *list);

/* pointer to the internal twobitencoding representation;
   it must be called after <gt_reads2twobit_encode> */
GtTwobitencoding *gt_reads2twobit_export_twobitencoding(
    const GtReads2Twobit *r2t);

/* pointer to the internal seppos array;
   it must be called after <gt_reads2twobit_encode> */
unsigned long *gt_reads2twobit_export_seppos(const GtReads2Twobit *r2t);

/* sort the sequences according to the specified comparator function <cmp> */
void gt_reads2twobit_sort(GtReads2Twobit *r2t, GtCompareWithData cmp,
    void *cmp_data);

/* write the libraries information to disk */
int gt_reads2twobit_write_libraries_table(const GtReads2Twobit *r2t,
    char *path, GtError *err);

#define GT_READS2TWOBIT_LIBSPEC_HELPMSG \
  "specify a list of input libraries (Fasta/FastQ); for single-end " \
  "libraries use the filename (which is not allowed to contain ':' " \
  "symbols); for paired-end libraries with reads interleaved (f,r,f,r,...) "\
  "in a single file use the notation <filename>:<insertlength>[,<stdev>] " \
  "(stdev may be omitted); for paired-end with reads in two files (f, r) " \
  "use the notation <file_f>:<file_r>:<insertlength>[,<stdev>]"

#endif
