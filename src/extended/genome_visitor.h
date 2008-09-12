/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef GENOME_VISITOR_H
#define GENOME_VISITOR_H

/* the ``genome visitor'' interface, a visitor for genome nodes */
typedef struct GenomeVisitorClass GenomeVisitorClass;
typedef struct GenomeVisitor GenomeVisitor;

#include "extended/comment.h"
#include "extended/genome_feature.h"
#include "extended/sequence_region.h"
#include "extended/sequence_node.h"

int  genome_visitor_visit_comment(GenomeVisitor*, GT_Comment*, GtError*);
int  genome_visitor_visit_genome_feature(GenomeVisitor*, GT_GenomeFeature*,
                                         GtError*);
int  genome_visitor_visit_sequence_region(GenomeVisitor*, GT_SequenceRegion*,
                                          GtError*);
int  genome_visitor_visit_sequence_node(GenomeVisitor*, GT_SequenceNode*,
                                        GtError*);
void genome_visitor_delete(GenomeVisitor *gv);

#endif
