/*
  Copyright (c) 2009-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef DUP_FEATURE_VISITOR_H
#define DUP_FEATURE_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtDupFeatureVisitor GtDupFeatureVisitor;

#include "extended/node_visitor.h"

const GtNodeVisitorClass* gt_dup_feature_visitor_class(void);
/* Duplicate internal feature nodes of type <source_type> as features with type
   <dest_type>. The duplicated feature does not inherit the children. */
GtNodeVisitor*            gt_dup_feature_visitor_new(const char *dest_type,
                                                     const char *source_type);

#endif
