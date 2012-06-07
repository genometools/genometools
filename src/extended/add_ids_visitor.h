/*
  Copyright (c) 2010-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef ADD_IDS_VISITOR_H
#define ADD_IDS_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtAddIDsVisitor GtAddIDsVisitor;

#include "extended/node_visitor.h"

const GtNodeVisitorClass* gt_add_ids_visitor_class(void);
GtNodeVisitor* gt_add_ids_visitor_new(bool ensure_sorting);
unsigned long  gt_add_ids_visitor_node_buffer_size(GtNodeVisitor*);
GtGenomeNode*  gt_add_ids_visitor_get_node(GtNodeVisitor*);
void           gt_add_ids_visitor_finalize(GtNodeVisitor*);

#endif
