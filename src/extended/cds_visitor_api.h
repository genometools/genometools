/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef CDS_VISITOR_API_H
#define CDS_VISITOR_API_H

/* <GtCDSVisitor> is a <GtNodeVisitor> that adds CDS features for the
   longest ORFs in a <GtFeatureNode>. */
typedef struct GtCDSVisitor GtCDSVisitor;

#include "extended/node_visitor_api.h"
#include "extended/region_mapping_api.h"

const GtNodeVisitorClass* gt_cds_visitor_class(void);

/* Create a new <GtCDSVisitor> with the given <region_mapping> for sequence,
   <minorflen> minimum ORF length, <source> as the source string.
   If <start_codon> is <true> a frame has to start with a start codon,
   otherwise a frame can start everywhere (i.e., at the first amino acid or
   after a stop codon). If <final_stop_codon> is <true> the last ORF must
   end with a stop codon, otherwise it can be ``open''.
   Takes ownership of <region_mapping>. */
GtNodeVisitor*            gt_cds_visitor_new(GtRegionMapping *region_mapping,
                                             unsigned int minorflen,
                                             GtStr *source, bool start_codon,
                                             bool final_stop_codon,
                                             bool generic_start_codons);

/* Sets <region_mapping> to be the region mapping specifying the sequence for
   <cds_visitor>. Does not take ownership of <region_mapping>. */
void                      gt_cds_visitor_set_region_mapping(GtCDSVisitor
                                                            *cds_visitor,
                                                            GtRegionMapping
                                                            *region_mapping);

#define                   gt_cds_visitor_cast(GV)\
                          gt_node_visitor_cast(gt_cds_visitor_class(), GV)

#endif
