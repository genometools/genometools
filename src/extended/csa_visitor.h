/*
  Copyright (c) 2006-2007, 2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007       Center for Bioinformatics, University of Hamburg

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

#ifndef CSA_VISITOR_H
#define CSA_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct CSAVisitor CSAVisitor;

#include "extended/node_visitor.h"

const GtNodeVisitorClass* gt_csa_visitor_class(void);
GtNodeVisitor*            gt_csa_visitor_new(unsigned long join_length);
unsigned long             gt_csa_visitor_node_buffer_size(GtNodeVisitor*);
GtGenomeNode*             gt_csa_visitor_get_node(GtNodeVisitor*);
void                      gt_csa_visitor_process_cluster(GtNodeVisitor*,
                                                         bool final_cluster);

#endif
