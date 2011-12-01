/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef ESA_VISITOR_REP_H
#define ESA_VISITOR_REP_H

#include <stdio.h>
#include "match/esa_visitor.h"

typedef void              (*GtESAVisitorFreeFunc)(GtESAVisitor*);
typedef GtESAVisitorInfo* (*GtESAVisitorInfoCreatorFunc)(GtESAVisitor*);
typedef void              (*GtESAVisitorInfoDestructorFunc)(GtESAVisitorInfo*,
                                                            GtESAVisitor*);
typedef int               (*GtESAVisitorLeafEdgeFunc)(GtESAVisitor*,
                                                      bool,
                                                      unsigned long,
                                                      unsigned long,
                                                      GtESAVisitorInfo*,
                                                      unsigned long,
                                                      GtError*);
typedef int               (*GtESAVisitorBranchingEdgeFunc)(GtESAVisitor*,
                                                           bool,
                                                           unsigned long,
                                                           unsigned long,
                                                           GtESAVisitorInfo*,
                                                           unsigned long,
                                                           unsigned long,
                                                           unsigned long,
                                                           GtESAVisitorInfo*,
                                                           GtError*);
typedef int               (*GtESAVisitorLCPIntervalFunc)(GtESAVisitor*,
                                                         unsigned long,
                                                         unsigned long,
                                                         unsigned long,
                                                         GtESAVisitorInfo*,
                                                         GtError*);

typedef struct GtESAVisitorMembers GtESAVisitorMembers;

struct GtESAVisitor {
  const GtESAVisitorClass *c_class;
  GtESAVisitorMembers *members;
};

const GtESAVisitorClass* gt_esa_visitor_class_new(size_t size,
                                                GtESAVisitorFreeFunc,
                                                GtESAVisitorLeafEdgeFunc,
                                                GtESAVisitorBranchingEdgeFunc,
                                                GtESAVisitorLCPIntervalFunc,
                                                GtESAVisitorInfoCreatorFunc,
                                                GtESAVisitorInfoDestructorFunc);
GtESAVisitor*            gt_esa_visitor_create(const GtESAVisitorClass*);

#endif
