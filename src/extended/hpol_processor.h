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

#ifndef HPOL_PROCESSOR_H
#define HPOL_PROCESSOR_H

#include "core/encseq.h"
#include "extended/seqpos_classifier.h"
#include "extended/aligned_segments_pile.h"
#include "core/logger.h"

/* Homopolymer processor; scans a sequence and searches homopolymers
 * (defined as stretches of at least <hmin> identical symbols) and
 * processes them. */

typedef struct GtHpolProcessor GtHpolProcessor;

GtHpolProcessor *gt_hpol_processor_new(GtEncseq *encseq, unsigned long hmin);
/* Enable the correction of length of the homopolymers in the segments
 * provided by <asp>. The <max_hlen_diff> parameter defines the maximal
 * difference in length from the reference for a correction to take place.
 * The <min_alt_consensus> parameter defines the minimal alternative
 * consensus among segments in homopolymer length, by which no correction
 * is done: a value of 0.5 means e.g. 50% consensus. To disable the
 * feature use a value over 1.0. */
void gt_hpol_processor_enable_segments_hlen_adjustment(GtHpolProcessor *hpp,
    GtAlignedSegmentsPile *asp, unsigned long max_hlen_diff,
    double min_alt_consensus);

/* Test: compare refregion of the segments with the reference sequence. */
void gt_hpol_processor_enable_aligned_segments_refregionscheck(
    GtHpolProcessor *hpp, GtAlignedSegmentsPile *asp);

/* Run the scanning over the sequence. If <logger> is not NULL, then verbose
 * information is output, including the distribution of homopolymers length.
 * Returns 0 on success, on error a negative number and <err> is set. */
int gt_hpol_processor_run(GtHpolProcessor *hpp, GtLogger *logger, GtError *err);

/* Restrict processing to homopolymers whose end position is classified
 * by <spc> as "inside" (e.g. inside CDS). */
void gt_hpol_processor_restrict_to_feature_type(GtHpolProcessor *hpp,
    GtSeqposClassifier *spc);

/* Output the segments sequence in FastQ format to <outfile> during run. */
void gt_hpol_processor_enable_segments_output(GtHpolProcessor *hpp,
    GtFile *outfile);

void gt_hpol_processor_delete(GtHpolProcessor *hpp);

#endif
