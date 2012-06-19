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

/* Homopolymer length corrector */

typedef struct GtHpolProcessor GtHpolProcessor;

GtHpolProcessor *gt_hpol_processor_new(GtEncseq *encseq, unsigned long hmin);

int gt_hpol_processor_run(GtHpolProcessor *hpp, GtLogger *logger, GtError *err);

void gt_hpol_processor_restrict_to_feature_type(GtHpolProcessor *hpp,
    GtSeqposClassifier *spc);

void gt_hpol_processor_enable_segments_hlen_adjustment(GtHpolProcessor *hpp,
    GtAlignedSegmentsPile *asp, unsigned long max_hlen_diff,
    double min_alt_consensus);

void gt_hpol_processor_enable_aligned_segments_refregionscheck(
    GtHpolProcessor *hpp, GtAlignedSegmentsPile *asp);

void gt_hpol_processor_enable_segments_output(GtHpolProcessor *hpp,
    GtFile *outfile);

void gt_hpol_processor_delete(GtHpolProcessor *hpp);

#endif
