/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef LTR_VISITOR_H
#define LTR_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtLTRVisitor GtLTRVisitor;

#include "extended/node_visitor.h"
#include "ltr/ltrelement.h"

const GtNodeVisitorClass* gt_ltr_visitor_class(void);
GtNodeVisitor*            gt_ltr_visitor_new(GtLTRElement *element);

#define gt_ltr_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltr_visitor_class(), GV)

#define gt_ltr_visitor_try_cast(GV)\
        gt_node_visitor_try_cast(gt_ltr_visitor_class(), GV)

#endif
