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

/* AlignedSegment contains an aligned <segment> (e.g. a read) to a <reference>.
 * The segment sequence (<seq>) and quality scores (<qual>) and the
 * aligned region of the reference sequence (this portion will be called here
 * <refregion>) are explicitely stored as NULL terminated editable strings. */
typedef struct GtAlignedSegment GtAlignedSegment;

/* create a new GtAlignedSegment instance fetching sequence and alignment
 * information from a SAM alignment. The <refregion> will contain
 * '?' at all positions which cannot be inferred from the data available
 * in the SAM alignment. If an unmapped SAM alignment is used,
 * the refregion will be NULL. The <sem> mapping will be used to
 * calculate the coordinates of refregion on the reference. */
GtAlignedSegment *gt_aligned_segment_new_from_sa(GtSamAlignment *sa,
    GtSamfileEncseqMapping *sem);

/* sequence of the segment */
char *gt_aligned_segment_seq(GtAlignedSegment *as);

/* quality scores of the segment */
char *gt_aligned_segment_qual(GtAlignedSegment *as);

/* sequence of <refregion> or NULL if segment is unmapped */
char *gt_aligned_segment_refregion(GtAlignedSegment *as);

/* make a copy of the original status of the segment sequence,
 * so that an edit script can be computed by comparing with the
 * edited sequence; this method may be called only once */
void gt_aligned_segment_enable_edit_tracking(GtAlignedSegment *as);

/* the original (aligned) status of the segment sequence;
 * if <gt_aligned_enable_edit_tracking> was not called, it returns NULL */
const char *gt_aligned_segment_orig_seq(GtAlignedSegment *as);

/* length of the unedited segment sequence */
unsigned long gt_aligned_segment_orig_seqlen(const GtAlignedSegment *as);

/* coordinate on the unedited ungapped segment sequence corresponding to
 * <refpos> on the reference sequence; returns the correct original segment
 * coordinates also in the case when the segment is reverse;
 * returns GT_UNDEF_ULONG if the coordinates are outside the refregion or
 * the segment is unmapped */
unsigned long gt_aligned_segment_orig_seqpos_for_refpos(
    const GtAlignedSegment *as, unsigned long refpos);

/* coordinates of <refregion> on the reference sequence
 * or GT_UNDEF_ULONG if unmapped */
unsigned long gt_aligned_segment_refregion_startpos(const GtAlignedSegment *as);
unsigned long gt_aligned_segment_refregion_endpos(const GtAlignedSegment *as);

/* finds the correct offset in the (eventually gapped) reference region
 * corresponding to a given reference sequence coordinate;
 * returns GT_UNDEF_ULONG if the coordinates are outside the region
 * or the segment is unmapped */
unsigned long gt_aligned_segment_offset_for_refpos(const GtAlignedSegment *as,
    unsigned long refpos);

/* total length of the alignment */
unsigned long gt_aligned_segment_length(const GtAlignedSegment *as);

/* remove gaps from <refregion>; assumes segment is not unmapped;
 * (after this operation segment and reference are not aligned anymore) */
void gt_aligned_segment_ungap_refregion(GtAlignedSegment *as);

/* remove gaps from the segment sequence and quality scores
 * (after this operation segment and reference are not aligned anymore) */
void gt_aligned_segment_ungap_seq_and_qual(GtAlignedSegment *as);

/* true if the alignment has insertions or deletions */
bool gt_aligned_segment_has_indels(const GtAlignedSegment *as);

/* true if the original alignment is on the opposite strand than the reference
 * sequence; the editable segment and refregion are however stored based the
 * reference strand sequence */
bool gt_aligned_segment_is_reverse(const GtAlignedSegment *as);

/* description line of the segment */
const char *gt_aligned_segment_description(const GtAlignedSegment *as);

/* change '?' symbols in refregion to real characters from an encseq;
 * (this only works correctly if refregion has either not been edited,
 * or at least the number of indels in it is balanced); it assumes the
 * segment is not unmapped; */
void gt_aligned_segment_assign_refregion_chars(GtAlignedSegment *as,
    GtEncseq *encseq);

/* show the current segment and refregion to file */
void gt_aligned_segment_show(GtAlignedSegment *as, GtFile *outfp);

/* mapping quality of the segment (GT_UNDEF_ULONG if unmapped) */
unsigned long gt_aligned_segment_mapping_quality(GtAlignedSegment *as);

/* getter/setter for edited bit of seq and refregion */
void gt_aligned_segment_seq_set_edited(GtAlignedSegment *as);
bool gt_aligned_segment_seq_edited(const GtAlignedSegment *as);
void gt_aligned_segment_refregion_set_edited(GtAlignedSegment *as);
bool gt_aligned_segment_refregion_edited(const GtAlignedSegment *as);

/* delete the aligned segment */
void gt_aligned_segment_delete(GtAlignedSegment *as);

#endif
