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

#ifndef SAMFILE_ENCSEQ_MAPPING_H
#define SAMFILE_ENCSEQ_MAPPING_H

#include "core/encseq.h"
#include "core/error_api.h"
#include "extended/samfile_iterator.h"

/* The <GtSamfileEncseqMapping> class represents a mapping of reference
   sequences of a SAM file to a <GtEncseq> containing the same sequences. */
typedef struct GtSamfileEncseqMapping GtSamfileEncseqMapping;

/* Returns a new <GtSamfileEncseqMapping> for <samfile_iterator> and <encseq>.
   The set of sequences encoded in the <encseq> and the set of references
   in the <samfile_iterator> must be identical (two sequences are hereby
   considered identical if the sequence identifiers are equal, i.e. the
   description up to the first space or newline).
   If the sequence sets are not identical, NULL is returned and <err> is set
   accordingly. */
GtSamfileEncseqMapping* gt_samfile_encseq_mapping_new(
                                            GtSamfileIterator *samfile_iterator,
                                            const GtEncseq *encseq,
                                            GtError *err);

/* Returns the coordinate in the <GtEncseq> corresponding to a given read with
   number <reference_num> in the SAM file and <reference_seqpos> */
unsigned long           gt_samfile_encseq_mapping_seqpos(
                                 GtSamfileEncseqMapping *samfile_encseq_mapping,
                                 int32_t reference_num,
                                 unsigned long reference_seqpos);

/* Deletes <samfile_encseq_mapping> and frees all associated space. */
void                    gt_samfile_encseq_mapping_delete(
                                GtSamfileEncseqMapping *samfile_encseq_mapping);

#endif
