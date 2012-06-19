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
 * Iterating over a sorted set of positions <position> in a reference sequence
 * determines a pile of segment sequences (e.g. reads), which are mapped
 * over <position>, according to a given mapping (SAM or BAM format).
 *
 * requirements/limitations:
 * - the SAM/BAM sorting order must be "coordinate"
 *   (e.g. output of samtools sort)
 */

typedef struct GtAlignedSegmentsPile GtAlignedSegmentsPile;

GtAlignedSegmentsPile *gt_aligned_segments_pile_new(GtSamfileIterator *sfi);

void gt_aligned_segments_pile_delete(GtAlignedSegmentsPile *asp);

/* move the pile over the next position, which must be
 * larger than the previous position, if the function
 * has been already called before */
void gt_aligned_segments_pile_move_over_position(
    GtAlignedSegmentsPile *asp, unsigned long position);

/* returns the pile as GtDlist of GtAlignedSegment objects,
 * which is sorted by end coordinate on the ref sequence */
const GtDlist *gt_aligned_segments_pile_get(const GtAlignedSegmentsPile *asp);

/* number of segments currently on the pile */
unsigned long gt_aligned_segments_pile_size(GtAlignedSegmentsPile *asp);

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

/* move away from current position (i.e. calling process_complete
 * on segments on the pile) and, if skip_remaining, call
 * process_unmapped / process_skipped on all remaining alignments
 * in the SAM iterator */
void gt_aligned_segments_pile_flush(GtAlignedSegmentsPile *asp,
    bool skip_remaining);

#endif
