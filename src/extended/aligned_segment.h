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

#ifndef ALIGNED_SEGMENT_H
#define ALIGNED_SEGMENT_H

#include "core/encseq.h"
#include "core/file.h"
#include "extended/sam_alignment.h"
#include "extended/samfile_encseq_mapping.h"

/* GtAlignedSegment contains a read ("segment") aligned to a
   region of a reference sequence. The segment sequence
   and quality scores are NULL terminated editable strings. */
typedef struct GtAlignedSegment GtAlignedSegment;

/* Creates a new GtAlignedSegment instance fetching sequence and alignment
   information from a SAM alignment. The <sem> mapping will be used to
   calculate the coordinates of refregion on the reference. */
GtAlignedSegment *gt_aligned_segment_new_from_sa(GtSamAlignment *sa,
    GtSamfileEncseqMapping *sem);

/* sequence of the segment */
char *gt_aligned_segment_seq(GtAlignedSegment *as);

/* quality scores of the segment */
char *gt_aligned_segment_qual(GtAlignedSegment *as);

/* sequence of the reference region */
char *gt_aligned_segment_refregion(GtAlignedSegment *as);

/* description of the segment */
const char *gt_aligned_segment_description(const GtAlignedSegment *as);

/* original length of the ungapped segment sequence */
unsigned long gt_aligned_segment_orig_seqlen(const GtAlignedSegment *as);

/* current length of the alignment (i.e. including gaps) */
unsigned long gt_aligned_segment_length(const GtAlignedSegment *as);

/* true if the original alignment is on the opposite strand than the reference
 * sequence; the editable segment and refregion are however stored based the
 * reference strand sequence */
bool gt_aligned_segment_is_reverse(const GtAlignedSegment *as);

/* true if the alignment has insertions or deletions */
bool gt_aligned_segment_has_indels(const GtAlignedSegment *as);

/* mapping quality of the segment */
unsigned long gt_aligned_segment_mapping_quality(GtAlignedSegment *as);

/* sets the edited bit of <as> */
void gt_aligned_segment_seq_set_edited(GtAlignedSegment *as);

/* gets the edited bit of <as> */
bool gt_aligned_segment_seq_edited(const GtAlignedSegment *as);

/* starting coordinate of the reference region on the reference sequence */
unsigned long gt_aligned_segment_refregion_startpos(const GtAlignedSegment *as);

/* ending coordinates of the reference region on the reference sequence */
unsigned long gt_aligned_segment_refregion_endpos(const GtAlignedSegment *as);

/* coordinates in the reference region of a given reference coordinate
   taking gaps into account */
unsigned long gt_aligned_segment_offset_for_refpos(const GtAlignedSegment *as,
    unsigned long refpos);

/* changes '?' symbols in refregion of <as> to real characters from <encseq>;
   useful for debug output */
void gt_aligned_segment_assign_refregion_chars(GtAlignedSegment *as,
    GtEncseq *encseq);

/* remove gaps from <refregion> for output */
void gt_aligned_segment_ungap_refregion(GtAlignedSegment *as);

/* remove gaps from the segment sequence and quality scores for output */
void gt_aligned_segment_ungap_seq_and_qual(GtAlignedSegment *as);

/* output the current segment and refregion to the file <outfp> */
void gt_aligned_segment_show(GtAlignedSegment *as, GtFile *outfp);

/* stores a copy of the original status of the segment sequence of <as>,
   allowing to compute an edit script by comparing with the edited sequence;
   (this is currently used by HOP -stats option, but can be done more
   efficiently) */
void gt_aligned_segment_enable_edit_tracking(GtAlignedSegment *as);

/* the original (aligned) status of the segment sequence;
 * needs <gt_aligned_enable_edit_tracking> */
const char *gt_aligned_segment_orig_seq(GtAlignedSegment *as);

/* coordinate in the unedited segment sequence corresponding to <refpos>
   on the reference sequence; takes gaps into account and works correctly
   also if the alignment is reverse; needs <gt_aligned_enable_edit_tracking> */
unsigned long gt_aligned_segment_orig_seqpos_for_refpos(
    const GtAlignedSegment *as, unsigned long refpos);

/* delete the aligned segment */
void gt_aligned_segment_delete(GtAlignedSegment *as);

#endif
