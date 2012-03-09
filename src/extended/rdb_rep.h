/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef RDB_REP_H
#define RDB_REP_H

#include <stdio.h>
#include "core/cstr_table_api.h"
#include "extended/rdb_api.h"

struct GtRDB {
  const GtRDBClass *c_class;
  GtRDBMembers *members;
};

struct GtRDBStmt {
  const GtRDBStmtClass *c_class;
};

typedef void (*GtRDBFreeFunc)(GtRDB*);
typedef GtRDBStmt* (*GtRDBPrepareFunc)(GtRDB*, const char*, unsigned long,
                                       GtError*);
typedef unsigned long (*GtRDBGetLastInsertIDFunc)(GtRDB*, const char*,
                                                  GtError*);
typedef int (*GtRDBAcceptVisitorFunc)(GtRDB*, GtRDBVisitor*, GtError*);
typedef GtCstrTable* (*GtRDBGetIndexesFunc)(GtRDB*, GtError*);
typedef GtCstrTable* (*GtRDBGetTablesFunc)(GtRDB*, GtError*);
typedef int (*GtRDBRecreateFunc)(GtRDB*, GtError*);

typedef int  (*GtRDBStmtResetFunc)(GtRDBStmt *stmt, GtError *err);
typedef int  (*GtRDBStmtBindIntFunc)(GtRDBStmt*, unsigned long, int, GtError*);
typedef int  (*GtRDBStmtBindUlongFunc)(GtRDBStmt*, unsigned long, unsigned long,
                                       GtError*);
typedef int  (*GtRDBStmtBindStringFunc)(GtRDBStmt*, unsigned long, const char*,
                                        GtError*);
typedef int  (*GtRDBStmtBindDoubleFunc)(GtRDBStmt*, unsigned long, double,
                                        GtError*);
typedef int  (*GtRDBStmtFetchFunc)(GtRDBStmt*, GtError*);
typedef void (*GtRDBStmtFreeFunc)(GtRDBStmt*);
typedef int  (*GtRDBStmtGetIntFunc)(GtRDBStmt*, unsigned long, int*, GtError*);
typedef int  (*GtRDBStmtGetUlongFunc)(GtRDBStmt*, unsigned long, unsigned long*,
                                     GtError*);
typedef int  (*GtRDBStmtGetStringFunc)(GtRDBStmt*, unsigned long, GtStr*,
                                      GtError*);
typedef int  (*GtRDBStmtGetDoubleFunc)(GtRDBStmt*, unsigned long, double*,
                                      GtError*);

const GtRDBClass* gt_rdb_class_new(size_t size,
                                   GtRDBFreeFunc free_func,
                                   GtRDBPrepareFunc prepare_func,
                                   GtRDBGetLastInsertIDFunc last_id_func,
                                   GtRDBAcceptVisitorFunc accept_func,
                                   GtRDBGetIndexesFunc get_indexes_func,
                                   GtRDBGetTablesFunc get_tables_func);

const GtRDBStmtClass* gt_rdb_stmt_class_new(size_t size,
                                       GtRDBStmtResetFunc reset_func,
                                       GtRDBStmtBindIntFunc bind_int_func,
                                       GtRDBStmtBindUlongFunc bind_ulong_func,
                                       GtRDBStmtBindStringFunc bind_string_func,
                                       GtRDBStmtBindDoubleFunc bind_double_func,
                                       GtRDBStmtFetchFunc fetch_func,
                                       GtRDBStmtGetIntFunc get_int_func,
                                       GtRDBStmtGetUlongFunc get_ulong_func,
                                       GtRDBStmtGetStringFunc get_string_func,
                                       GtRDBStmtGetDoubleFunc get_double_func,
                                       GtRDBStmtFreeFunc free_func);

GtRDB*     gt_rdb_create(const GtRDBClass*);
GtRDBStmt* gt_rdb_stmt_create(const GtRDBStmtClass*);
void*      gt_rdb_cast(const GtRDBClass*, GtRDB*);
void*      gt_rdb_stmt_cast(const GtRDBStmtClass*, GtRDBStmt*);

#endif
