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

#ifndef ALIGNED_SEGMENTS_PILE_H
#define ALIGNED_SEGMENTS_PILE_H

#include "core/dlist_api.h"
#include "extended/samfile_iterator.h"
#include "extended/aligned_segment.h"

/*
 * A GtAlignedSegmentsPile represents a pile of GtAlignedSegment objects,
 * segment sequences (e.g. reads), which are mapped over <position>,
 * according to a given mapping file (SAM or BAM format).
 *
 * The underlying mapping file is accessed sequentially by a GtSamfileIterator
 * (which is specifed by calling the <gt_aligned_segments_pile_new> method).
 * As access is sequential, the positions over which a pile is constructed
 * must be visited in sorted order (in other words: piles can only be
 * constructed over positions which are after previously visited positions).
 *
 * A requirement/limitation is that the SAM/BAM sorting order must be
 * "coordinate" (e.g. output of samtools sort).
 *
 */

typedef struct GtAlignedSegmentsPile GtAlignedSegmentsPile;

GtAlignedSegmentsPile *gt_aligned_segments_pile_new(GtSamfileIterator *sfi);

/* moves the pile over <position>, which must be larger than any previous
 * position (if the function has been already called before) */
void gt_aligned_segments_pile_move_over_position(
    GtAlignedSegmentsPile *asp, unsigned long position);

/* a double linked list (GtDlist) of GtAlignedSegment objects,
 * sorted by end coordinate on the reference sequence */
const GtDlist *gt_aligned_segments_pile_get(const GtAlignedSegmentsPile *asp);

/* number of segments currently on the pile */
unsigned long gt_aligned_segments_pile_size(GtAlignedSegmentsPile *asp);

/* enable edit tracking feature of the aligned segments */
void gt_aligned_segments_pile_enable_edit_tracking(GtAlignedSegmentsPile *asp);

/* processing function type used by the hook functions */
typedef void (*GtAlignedSegmentsPileProcessFunc)(GtAlignedSegment *as,
    void *data);

/* register an hook, which is called each time an alignment rightmost
 * position has been reached and the alignment will thus be removed
 * from the pile */
void gt_aligned_segments_pile_register_process_complete(
    GtAlignedSegmentsPile *asp,
    GtAlignedSegmentsPileProcessFunc process_complete, void *pc_data);

/* register an hook, which is called each time the SamfileIterator
 * visits a read which lies completely between the previous
 * and the current <position> of the pile */
void gt_aligned_segments_pile_register_process_skipped(
    GtAlignedSegmentsPile *asp,
    GtAlignedSegmentsPileProcessFunc process_skipped, void *ps_data);

/* register an hook, which is called each time the SamfileIterator
 * visits an unmapped read */
void gt_aligned_segments_pile_register_process_unmapped(
    GtAlignedSegmentsPile *asp,
    GtAlignedSegmentsPileProcessFunc process_unmapped, void *pu_data);

/* prepare the pile for deletion, by moving away from current position
 * (process_complete is called on segments on the pile);
 * furthmore if skip_remaining, call process_unmapped / process_skipped
 * is called for all remaining alignments visited by the GtSamIterator */
void gt_aligned_segments_pile_flush(GtAlignedSegmentsPile *asp,
    bool skip_remaining);

/* disable deletion of aligned sequence segments after processing;
 * if this method is called, the user must manually delete each segment
 * (thus the processing functions hooks are not allowed to be NULL) */
void gt_aligned_segments_pile_disable_segment_deletion(
    GtAlignedSegmentsPile *asp);

/* flush and delete the pile */
void gt_aligned_segments_pile_delete(GtAlignedSegmentsPile *asp);

#endif
