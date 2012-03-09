/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef RDB_VISITOR_REP_H
#define RDB_VISITOR_REP_H

#include "core/error_api.h"
#include "extended/rdb_visitor_api.h"

typedef void (*GtRDBVisitorFreeFunc)(GtRDBVisitor*);
typedef int  (*GtRDBVisitorSqliteFunc)(GtRDBVisitor*, GtRDBSqlite*, GtError*);
typedef int  (*GtRDBVisitorMySQLFunc)(GtRDBVisitor*, GtRDBMySQL*, GtError*);

typedef struct GtRDBVisitorClass GtRDBVisitorClass;
typedef struct GtRDBVisitorMembers GtRDBVisitorMembers;

struct GtRDBVisitor {
  const GtRDBVisitorClass *c_class;
  GtRDBVisitorMembers *members;
};

const GtRDBVisitorClass* gt_rdb_visitor_class_new(size_t size,
                                                  GtRDBVisitorFreeFunc,
                                                  GtRDBVisitorSqliteFunc,
                                                  GtRDBVisitorMySQLFunc);
GtRDBVisitor*       gt_rdb_visitor_create(const GtRDBVisitorClass*);
void*               gt_rdb_visitor_cast(const GtRDBVisitorClass*,
                                         GtRDBVisitor*);

#endif
