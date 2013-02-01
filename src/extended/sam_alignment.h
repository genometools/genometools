/*
  Copyright (c) 2011      Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c)      2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef SAM_ALIGNMENT_H
#define SAM_ALIGNMENT_H

#ifndef S_SPLINT_S
#include <stdbool.h>
#endif

#include "core/alphabet_api.h"
#include "core/str_api.h"
#include "core/types_api.h"

typedef struct GtSamAlignment GtSamAlignment;

GtSamAlignment* gt_sam_alignment_new(GtAlphabet *alphabet);

const char*     gt_sam_alignment_identifier(GtSamAlignment *sam_alignment);

GtSamAlignment* gt_sam_alignment_clone(GtSamAlignment *sam_alignment);

/* Returns the number of the reference sequence this alignment corresponds to,
   mapping of numbers to names is done in the sam/bam-header. Access through
   GtSamfileIterator */
int32_t         gt_sam_alignment_ref_num(GtSamAlignment *sam_alignment);

/* Returns the starting position of the alignment in the reference sequence */
unsigned long   gt_sam_alignment_pos(GtSamAlignment *sam_alignment);

/* Returns the ending position of the alignment in the reference sequence */
unsigned long gt_sam_alignment_rightmost_pos(GtSamAlignment *sam_alignment);

/* Returns length of read, not length of the alignment */
unsigned long   gt_sam_alignment_read_length(GtSamAlignment *sam_alignment);

/* Returns mapping quality value */
unsigned long gt_sam_alignment_mapping_quality(GtSamAlignment *sam_alignment);

/* Returns encoded read sequence from <sam_alignment>. */
const GtUchar*  gt_sam_alignment_sequence(GtSamAlignment *sam_alignment);

/* lower level version of <gt_sam_alignment_sequence> using an external
 * buffer to save the sequence information */
void gt_sam_alignment_sequence_external_buffer(GtSamAlignment *sam_alignment,
    GtUchar **seq_buffer, unsigned long *bufsize);

/* Returns string of qualities in ASCII format as in Sanger FASTQ for the
   read sequence from <sam_alignment>.
   The length is the same as the length of the read sequence. */
const GtUchar*  gt_sam_alignment_qualitystring(GtSamAlignment *sam_alignment);

/* Returns the number of CIGAR operations in <sam_alignment>. */
uint16_t        gt_sam_alignment_cigar_length(GtSamAlignment *sam_alignment);

/* lower level version of <gt_sam_alignment_qualitystring> using an external
 * buffer to save the sequence information */
void gt_sam_alignment_qualitystring_external_buffer(
    GtSamAlignment *sam_alignment, GtUchar **qual_buffer,
    unsigned long *bufsize);

/* Returns the length of CIGAR operation <i> in <sam_alignment>. */
uint32_t        gt_sam_alignment_cigar_i_length(GtSamAlignment *sam_alignment,
                                                uint16_t i);

/* Returns the type of CIGAR operation <i> in <sam_alignment>.
   Type is one of [MIDNSHP=X] (see sam/bam format documentation for details) */
unsigned char   gt_sam_alignment_cigar_i_operation(
                                                  GtSamAlignment *sam_alignment,
                                                  uint16_t i);
/* For explanation of the flag and how to interpret see samtools
   documentation. */
uint32_t        gt_sam_alignment_flag(GtSamAlignment *sam_alignment);

/* Checks the flag and returns true if bit is set in flag of <sam_alignment>.
   See sam/bam fileformat documentation for explanation of meaning of bits. */
bool            gt_sam_alignment_is_paired(GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_is_proper_paired(
                                                 GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_is_unmapped(GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_mate_is_unmapped(
                                                 GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_is_reverse(GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_mate_is_reverse(GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_is_read1(GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_is_read2(GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_is_secondary(GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_has_qc_failure(GtSamAlignment *sam_alignment);
bool            gt_sam_alignment_is_optical_pcr_duplicate(
                                                 GtSamAlignment *sam_alignment);

void            gt_sam_alignment_delete(GtSamAlignment *sam_alignment);
#endif
