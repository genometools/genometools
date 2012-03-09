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

#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/rdb_api.h"
#include "extended/rdb_rep.h"

struct GtRDBClass {
  size_t size;
  GtRDBFreeFunc free_func;
  GtRDBPrepareFunc prepare_func;
  GtRDBGetLastInsertIDFunc last_id_func;
  GtRDBAcceptVisitorFunc accept_func;
  GtRDBGetIndexesFunc get_indexes_func;
  GtRDBGetTablesFunc get_tables_func;
  GtRDBRecreateFunc recreate_func;
  GtRDBMembers *members;
};

struct GtRDBStmtClass {
  size_t size;
  GtRDBStmtResetFunc reset_func;
  GtRDBStmtBindIntFunc bind_int_func;
  GtRDBStmtBindUlongFunc bind_ulong_func;
  GtRDBStmtBindStringFunc bind_string_func;
  GtRDBStmtBindDoubleFunc bind_double_func;
  GtRDBStmtFetchFunc fetch_func;
  GtRDBStmtGetIntFunc get_int_func;
  GtRDBStmtGetUlongFunc get_ulong_func;
  GtRDBStmtGetStringFunc get_string_func;
  GtRDBStmtGetDoubleFunc get_double_func;
  GtRDBStmtFreeFunc free_func;
};

struct GtRDBMembers {
  unsigned int reference_count;
};

const GtRDBClass* gt_rdb_class_new(size_t size,
                                   GtRDBFreeFunc free_func,
                                   GtRDBPrepareFunc prepare_func,
                                   GtRDBGetLastInsertIDFunc last_id_func,
                                   GtRDBAcceptVisitorFunc accept_func,
                                   GtRDBGetIndexesFunc get_indexes_func,
                                   GtRDBGetTablesFunc get_tables_func)
{
  GtRDBClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free_func = free_func;
  c_class->prepare_func = prepare_func;
  c_class->accept_func = accept_func;
  c_class->last_id_func = last_id_func;
  c_class->get_indexes_func = get_indexes_func;
  c_class->get_tables_func = get_tables_func;
  return c_class;
}

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
                                       GtRDBStmtFreeFunc free_func)
{
  GtRDBStmtClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->reset_func = reset_func;
  c_class->bind_int_func = bind_int_func;
  c_class->bind_ulong_func = bind_ulong_func;
  c_class->bind_string_func = bind_string_func;
  c_class->bind_double_func = bind_double_func;
  c_class->fetch_func = fetch_func;
  c_class->get_int_func = get_int_func;
  c_class->get_ulong_func = get_ulong_func;
  c_class->get_string_func = get_string_func;
  c_class->get_double_func = get_double_func;
  c_class->free_func = free_func;
  return c_class;
}

GtRDB* gt_rdb_create(const GtRDBClass *dbc)
{
  GtRDB *ns;
  gt_assert(dbc && dbc->size);
  ns = gt_calloc(1, dbc->size);
  ns->c_class = dbc;
  ns->members = gt_calloc(1, sizeof (GtRDBMembers));
  return ns;
}

GtRDB* gt_rdb_ref(GtRDB *ns)
{
  gt_assert(ns);
  ns->members->reference_count++;
  return ns;
}

void gt_rdb_delete(GtRDB *ns)
{
  if (!ns) return;
  if (ns->members->reference_count) {
    ns->members->reference_count--;
    return;
  }
  gt_assert(ns->c_class);
  if (ns->c_class->free_func) ns->c_class->free_func(ns);
  gt_free(ns->members);
  gt_free(ns);
}

void* gt_rdb_cast(GT_UNUSED const GtRDBClass *dbc, GtRDB *db)
{
  gt_assert(dbc && db && db->c_class == dbc);
  return db;
}

int gt_rdb_accept(GtRDB *db, GtRDBVisitor *v, GtError *err)
{
  gt_assert(db && db->c_class);
  if (db->c_class->accept_func)
    return db->c_class->accept_func(db, v, err);
  return 0;
}

unsigned long gt_rdb_last_inserted_id(GtRDB *db, const char *table, GtError *e)
{
  gt_assert(db && db->c_class);
  if (db->c_class->last_id_func)
    return db->c_class->last_id_func(db, table, e);
  return GT_UNDEF_ULONG;
}

GtCstrTable* gt_rdb_get_indexes(GtRDB *db, GtError *err)
{
  gt_assert(db && db->c_class);
  if (db->c_class->get_indexes_func)
    return db->c_class->get_indexes_func(db, err);
  return NULL;
}

GtCstrTable* gt_rdb_get_tables(GtRDB *db, GtError *err)
{
  gt_assert(db && db->c_class);
  if (db->c_class->get_tables_func)
    return db->c_class->get_tables_func(db, err);
  return NULL;
}

GtRDBStmt* gt_rdb_stmt_create(const GtRDBStmtClass *stc)
{
  GtRDBStmt *st;
  gt_assert(stc && stc->size);
  st = gt_calloc(1, stc->size);
  st->c_class = stc;
  return st;
}

void* gt_rdb_stmt_cast(GT_UNUSED const GtRDBStmtClass *stmtc, GtRDBStmt *stmt)
{
  gt_assert(stmtc && stmt && stmt->c_class == stmtc);
  return stmt;
}

GtRDBStmt* gt_rdb_prepare(GtRDB *db, const char *query,
                          unsigned long num_params, GtError *err)
{
  gt_assert(db && db->c_class);
  if (db->c_class->prepare_func)
    return db->c_class->prepare_func(db, query, num_params, err);
  return 0;
}

int gt_rdb_stmt_reset(GtRDBStmt *stmt, GtError *err)
{
  gt_assert(stmt && stmt->c_class);
  if (stmt->c_class->reset_func)
    return stmt->c_class->reset_func(stmt, err);
  return 0;
}

int gt_rdb_stmt_bind_int(GtRDBStmt *stmt, unsigned long param_no,
                         int val, GtError *err)
{
  gt_assert(stmt && stmt->c_class);
  if (stmt->c_class->bind_int_func)
    return stmt->c_class->bind_int_func(stmt, param_no, val, err);
  return 0;
}

int gt_rdb_stmt_bind_ulong(GtRDBStmt *stmt, unsigned long param_no,
                           unsigned long val, GtError *err)
{
  gt_assert(stmt && stmt->c_class);
  if (stmt->c_class->bind_ulong_func)
    return stmt->c_class->bind_ulong_func(stmt, param_no, val, err);
  return 0;
}

int gt_rdb_stmt_bind_string(GtRDBStmt *stmt, unsigned long param_no,
                            const char *val, GtError *err)
{
  gt_assert(stmt && stmt->c_class);
  if (stmt->c_class->bind_string_func)
    return stmt->c_class->bind_string_func(stmt, param_no, val, err);
  return 0;
}

int gt_rdb_stmt_bind_double(GtRDBStmt *stmt, unsigned long param_no,
                            double val, GtError *err)
{
  gt_assert(stmt && stmt->c_class);
  if (stmt->c_class->bind_double_func)
    return stmt->c_class->bind_double_func(stmt, param_no, val, err);
  return 0;
}

int gt_rdb_stmt_exec(GtRDBStmt *stmt, GtError *err)
{
  gt_assert(stmt && stmt->c_class);
  if (stmt->c_class->fetch_func)
    return stmt->c_class->fetch_func(stmt, err);
  return 0;
}

void gt_rdb_stmt_delete(GtRDBStmt *stmt)
{
  if (!stmt) return;
  gt_assert(stmt->c_class);
  if (stmt->c_class->free_func)
    stmt->c_class->free_func(stmt);
  gt_free(stmt);
}

int gt_rdb_stmt_get_ulong(GtRDBStmt *stmt, unsigned long field_no,
                          unsigned long *result, GtError *err)
{
  gt_assert(stmt && stmt->c_class && result);
  if (stmt->c_class->get_ulong_func)
    return stmt->c_class->get_ulong_func(stmt, field_no, result, err);
  return 0;
}

int gt_rdb_stmt_get_int(GtRDBStmt *stmt, unsigned long field_no,
                        int *result, GtError *err)
{
  gt_assert(stmt && stmt->c_class && result);
  if (stmt->c_class->get_int_func)
    return stmt->c_class->get_int_func(stmt, field_no, result, err);
  return 0;
}

int gt_rdb_stmt_get_string(GtRDBStmt *stmt, unsigned long field_no,
                           GtStr *result, GtError *err)
{
  gt_assert(stmt && stmt->c_class && result);
  if (stmt->c_class->get_string_func)
    return stmt->c_class->get_string_func(stmt, field_no, result, err);
  return 0;
}

int gt_rdb_stmt_get_double(GtRDBStmt *stmt, unsigned long field_no,
                           double *result, GtError *err)
{
  gt_assert(stmt && stmt->c_class && result);
  if (stmt->c_class->get_double_func)
    return stmt->c_class->get_double_func(stmt, field_no, result, err);
  return 0;
}
