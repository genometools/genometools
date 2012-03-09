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

#include <stdlib.h>
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/rdb_visitor_rep.h"

/* the <GtRDBVisitor> interface */
struct GtRDBVisitorClass {
  size_t size;
  GtRDBVisitorFreeFunc free;
  GtRDBVisitorSqliteFunc sqlite_func;
  GtRDBVisitorMySQLFunc mysql_func;
};

const GtRDBVisitorClass*
gt_rdb_visitor_class_new(size_t size,
                          GtRDBVisitorFreeFunc free,
                          GtRDBVisitorSqliteFunc sqlite_func,
                          GtRDBVisitorMySQLFunc mysql_func)
{
  GtRDBVisitorClass *c_class;
  gt_assert(size);
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->sqlite_func = sqlite_func;
  c_class->mysql_func = mysql_func;
  return c_class;
}

GtRDBVisitor* gt_rdb_visitor_create(const GtRDBVisitorClass *rdbvc)
{
  GtRDBVisitor *rdbv;
  gt_assert(rdbvc && rdbvc->size);
  rdbv = gt_calloc(1, rdbvc->size);
  rdbv->c_class = rdbvc;
  return rdbv;
}

void* gt_rdb_visitor_cast(GT_UNUSED const GtRDBVisitorClass *rdbvc,
                          GtRDBVisitor *rdbv)
{
  gt_assert(rdbvc && rdbv && rdbv->c_class == rdbvc);
  return rdbv;
}

int gt_rdb_visitor_visit_sqlite(GtRDBVisitor *rdbv, GtRDBSqlite *rdbs,
                                GtError *err)
{
  gt_error_check(err);
  gt_assert(rdbv && rdbs && rdbv->c_class);
  if (rdbv->c_class->sqlite_func)
    return rdbv->c_class->sqlite_func(rdbv, rdbs, err);
  return 0;
}

int gt_rdb_visitor_visit_mysql(GtRDBVisitor *rdbv, GtRDBMySQL *rdbm,
                               GtError *err)
{
  gt_error_check(err);
  gt_assert(rdbv && rdbm && rdbv->c_class);
  if (rdbv->c_class->mysql_func)
    return rdbv->c_class->mysql_func(rdbv, rdbm, err);
  return 0;
}

void gt_rdb_visitor_delete(GtRDBVisitor *rdbv)
{
  if (!rdbv) return;
  gt_assert(rdbv->c_class);
  if (rdbv->c_class->free)
    rdbv->c_class->free(rdbv);
  gt_free(rdbv);
}
