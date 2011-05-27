/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef CDS_VISITOR_H
#define CDS_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtCDSVisitor GtCDSVisitor;

#include "extended/node_visitor.h"
#include "extended/region_mapping_api.h"

const GtNodeVisitorClass* gt_cds_visitor_class(void);
/* Takes ownership of <region_mapping>. */
GtNodeVisitor*            gt_cds_visitor_new(GtRegionMapping *region_mapping,
                                             unsigned int minorflen,
                                             GtStr *source, bool start_codon,
                                             bool final_stop_codon,
                                             bool generic_start_codons);
/* Does not take ownership of <region_mapping>. */
void                      gt_cds_visitor_set_region_mapping(GtCDSVisitor
                                                            *cds_visitor,
                                                            GtRegionMapping
                                                            *region_mapping);

#define                   gt_cds_visitor_cast(GV)\
                          gt_node_visitor_cast(gt_cds_visitor_class(), GV)

#endif
