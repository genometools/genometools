/*
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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
#ifndef WTREE_REP_H
#define WTREE_REP_H

#include <stdlib.h>

#include "extended/wtree.h"

typedef struct GtWtreeClass GtWtreeClass;
typedef struct GtWtreeMembers GtWtreeMembers;

typedef GtWtreeSymbol (*GtWtreeAccessFunc)(GtWtree*, GtUword);
typedef GtUword (*GtWtreeRankFunc)(GtWtree*, GtUword,
                                         GtWtreeSymbol);
typedef GtUword (*GtWtreeSelectFunc)(GtWtree*, GtUword,
                                          GtWtreeSymbol);
typedef void (*GtWtreeDeleteFunc)(GtWtree*);

struct GtWtree {
  const GtWtreeClass *c_class;
  GtWtreeMembers     *members;
};

struct GtWtreeClass {
  GtWtreeAccessFunc access_func;
  GtWtreeDeleteFunc delete_func;
  GtWtreeRankFunc   rank_func;
  GtWtreeSelectFunc select_func;
  size_t            size;
};

struct GtWtreeMembers {
  GtUword length,
                num_of_symbols,
                refcount;
};

const GtWtreeClass* gt_wtree_class_new(size_t size,
                                       GtWtreeAccessFunc,
                                       GtWtreeRankFunc,
                                       GtWtreeSelectFunc,
                                       GtWtreeDeleteFunc);

GtWtree*            gt_wtree_create(const GtWtreeClass*);

void*               gt_wtree_cast(const GtWtreeClass*,
                                  GtWtree*);

#endif
