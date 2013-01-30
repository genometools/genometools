/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef CSTR_ITERATOR_REP_H
#define CSTR_ITERATOR_REP_H

#include "extended/cstr_iterator.h"

typedef struct GtCstrIteratorClass GtCstrIteratorClass;
typedef struct GtCstrIteratorMembers GtCstrIteratorMembers;

typedef int    (*GtCstrIteratorNextFunc)(GtCstrIterator*,
                                         const char**,
                                         GtError*);
typedef int    (*GtCstrIteratorResetFunc)(GtCstrIterator*, GtError*);
typedef void   (*GtCstrIteratorDeleteFunc)(GtCstrIterator*);

struct GtCstrIterator {
  const GtCstrIteratorClass *c_class;
  GtCstrIteratorMembers *members;
};

struct GtCstrIteratorClass {
  size_t size;
  GtCstrIteratorNextFunc next_func;
  GtCstrIteratorResetFunc reset_func;
  GtCstrIteratorDeleteFunc delete_func;
};

GtCstrIteratorClass* gt_cstr_iterator_class_new(size_t size,
                                                GtCstrIteratorNextFunc,
                                                GtCstrIteratorResetFunc,
                                                GtCstrIteratorDeleteFunc);

GtCstrIterator*      gt_cstr_iterator_create(const GtCstrIteratorClass*);

void*                gt_cstr_iterator_cast(const GtCstrIteratorClass*,
                                           GtCstrIterator*);
#endif
