/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtext/comment.h"
#include "libgtext/genome_feature.h"
#include "libgtext/sequence_region.h"

int  genome_visitor_visit_comment(GenomeVisitor*, Comment*, Error*);
int  genome_visitor_visit_genome_feature(GenomeVisitor*, GenomeFeature*,
                                         Error*);
int  genome_visitor_visit_sequence_region(GenomeVisitor*, SequenceRegion*,
                                          Error*);
void genome_visitor_delete(GenomeVisitor *gv);

#endif
