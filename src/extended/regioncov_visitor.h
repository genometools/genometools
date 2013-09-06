/*
  Copyright (c) 2007, 2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007       Center for Bioinformatics, University of Hamburg

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

#ifndef REGIONCOV_VISITOR_H
#define REGIONCOV_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtRegionCovVisitor GtRegionCovVisitor;

#include "extended/node_visitor.h"

const GtNodeVisitorClass* gt_regioncov_visitor_class(void);
GtNodeVisitor*            gt_regioncov_visitor_new(GtUword
                                                   max_feature_dist);
void                      gt_regioncov_visitor_show_coverage(GtNodeVisitor*);

#endif
