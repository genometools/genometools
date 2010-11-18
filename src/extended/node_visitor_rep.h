/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef NODE_VISITOR_REP_H
#define NODE_VISITOR_REP_H

#include <stdio.h>
#include "extended/node_visitor.h"

typedef void (*GtNodeVisitorFreeFunc)(GtNodeVisitor*);
typedef int  (*GtNodeVisitorCommentNodeFunc)(GtNodeVisitor*, GtCommentNode*,
                                             GtError*);
typedef int  (*GtNodeVisitorFeatureNodeFunc)(GtNodeVisitor*, GtFeatureNode*,
                                             GtError*);
typedef int  (*GtNodeVisitorRegionNodeFunc)(GtNodeVisitor*, GtRegionNode*,
                                            GtError*);
typedef int  (*GtNodeVisitorSequenceNodeFunc)(GtNodeVisitor*, GtSequenceNode*,
                                              GtError*);
typedef int  (*GtNodeVisitorEOFNodeFunc)(GtNodeVisitor*, GtEOFNode*, GtError*);

typedef struct GtNodeVisitorMembers GtNodeVisitorMembers;

struct GtNodeVisitor {
  const GtNodeVisitorClass *c_class;
  GtNodeVisitorMembers *members;
};

const
GtNodeVisitorClass* gt_node_visitor_class_new(size_t size,
                                              GtNodeVisitorFreeFunc,
                                              GtNodeVisitorCommentNodeFunc,
                                              GtNodeVisitorFeatureNodeFunc,
                                              GtNodeVisitorRegionNodeFunc,
                                              GtNodeVisitorSequenceNodeFunc,
                                              GtNodeVisitorEOFNodeFunc);
GtNodeVisitor*      gt_node_visitor_create(const GtNodeVisitorClass*);

#endif
