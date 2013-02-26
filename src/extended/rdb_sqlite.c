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

#ifdef HAVE_SQLITE

#include <sqlite3.h>
#include <string.h>
#include "core/array_api.h"
#include "core/class_alloc_lock.h"
#include "core/unused_api.h"
#include "extended/rdb_rep.h"
#include "extended/rdb_sqlite_api.h"

#define GT_SQLITE_ERRMSG  "SQLite error code %d: %s"

struct GtRDBSqlite {
  const GtRDB parent_instance;
  sqlite3 *db;
};

struct GtRDBStmtSqlite {
  const GtRDBStmt parent_instance;
  sqlite3_stmt *stmt;
  sqlite3 *db;
  unsigned long num_params;
};

const GtRDBClass*     gt_rdb_sqlite_class(void);
const GtRDBStmtClass* gt_rdb_stmt_sqlite_class(void);

#define gt_rdb_sqlite_cast(RDB)\
        gt_rdb_cast(gt_rdb_sqlite_class(), RDB)

#define gt_rdb_stmt_sqlite_cast(RDB)\
        gt_rdb_stmt_cast(gt_rdb_stmt_sqlite_class(), RDB)

GtRDB* gt_rdb_sqlite_new(const char *dbpath, GtError *err)
{
  GtRDBSqlite *rdbs;
  int retval, had_err = 0;
  sqlite3 *db = NULL;
  GtRDB *rdb = NULL;
  gt_assert(dbpath);
  gt_error_check(err);
  retval = sqlite3_open(dbpath, &db);
  if (retval != SQLITE_OK) {
    gt_error_set(err, "cannot open database: %s", sqlite3_errmsg(db));
    sqlite3_close(db);
    had_err = -1;
  }
  if (!had_err) {
    rdb = gt_rdb_create(gt_rdb_sqlite_class());
    rdbs = gt_rdb_sqlite_cast(rdb);
    rdbs->db = db;
    gt_assert(rdbs->db);
  }
  return (had_err ? NULL : rdb);
}

static unsigned long gt_rdb_sqlite_last_inserted_id(GtRDB *rdb,
                                                   GT_UNUSED const char *table,
                                                   GT_UNUSED GtError *err)
{
  GtRDBSqlite *rdbs;
  gt_assert(rdb);
  gt_error_check(err);
  rdbs = gt_rdb_sqlite_cast(rdb);
  return (unsigned long) sqlite3_last_insert_rowid(rdbs->db);
}

static int gt_rdb_sqlite_accept(GtRDB *rdb, GtRDBVisitor *v, GtError *err)
{
  GtRDBSqlite *rdbs;
  gt_assert(rdb && v);
  gt_error_check(err);
  rdbs = gt_rdb_sqlite_cast(rdb);
  return gt_rdb_visitor_visit_sqlite(v, rdbs, err);
}

static void gt_rdb_sqlite_delete(GtRDB *rdb)
{
  GtRDBSqlite *rdbs;
  if (!rdb) return;
  rdbs = gt_rdb_sqlite_cast(rdb);
  if (rdbs->db)
    sqlite3_close(rdbs->db);
}

static GtRDBStmt* gt_rdb_sqlite_prepare(GtRDB *rdb, const char *query,
                                        unsigned long num_params, GtError *err)
{
  GtRDBStmt *st = NULL;
  GtRDBStmtSqlite *sts = NULL;
  GtRDBSqlite *rdbs;
  int retval, had_err = 0;
  sqlite3_stmt *tmp = NULL;
  gt_assert(rdb && query);
  gt_error_check(err);
  rdbs = gt_rdb_sqlite_cast(rdb);
  retval = sqlite3_prepare_v2(rdbs->db, query, -1, &tmp, NULL);
  if (retval != SQLITE_OK) {
    gt_error_set(err, GT_SQLITE_ERRMSG, retval, sqlite3_errmsg(rdbs->db));
    had_err = -1;
  }
  if (!had_err) {
    st = gt_rdb_stmt_create(gt_rdb_stmt_sqlite_class());
    sts = gt_rdb_stmt_sqlite_cast(st);
    sts->num_params = num_params;
    sts->stmt = tmp;
    sts->db = rdbs->db;
  }
  return st;
}

static int gt_rdb_stmt_sqlite_reset(GtRDBStmt *st, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int rval, had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  rval = sqlite3_reset(sts->stmt);
  if (rval != SQLITE_OK) {
    gt_error_set(err, GT_SQLITE_ERRMSG, rval, sqlite3_errmsg(sts->db));
    had_err = -1;
  }
  if (!had_err) {
    rval = sqlite3_clear_bindings(sts->stmt);
    if (rval != SQLITE_OK) {
      gt_error_set(err, GT_SQLITE_ERRMSG, rval, sqlite3_errmsg(sts->db));
      had_err = -1;
    }
  }
  return had_err;
}

static void gt_rdb_stmt_sqlite_delete(GtRDBStmt *st)
{
  GtRDBStmtSqlite *sts;
  if (!st) return;
  sts = gt_rdb_stmt_sqlite_cast(st);
  if (sts->stmt)
    sqlite3_finalize(sts->stmt);
}

static int gt_rdb_stmt_sqlite_bind_int(GtRDBStmt *st, unsigned long param_no,
                                       int val, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int rval, had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  gt_assert(param_no < sts->num_params);
  rval = sqlite3_bind_int(sts->stmt, param_no + 1, val);
  if (rval != SQLITE_OK) {
    gt_error_set(err, GT_SQLITE_ERRMSG, rval, sqlite3_errmsg(sts->db));
    had_err = -1;
  }
  return had_err;
}

static int gt_rdb_stmt_sqlite_bind_ulong(GtRDBStmt *st, unsigned long param_no,
                                         unsigned long val, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int rval, had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  gt_assert(param_no < sts->num_params);
  rval = sqlite3_bind_int(sts->stmt, param_no + 1, (int) val);
  if (rval != SQLITE_OK) {
    gt_error_set(err, GT_SQLITE_ERRMSG, rval, sqlite3_errmsg(sts->db));
    had_err = -1;
  }
  return had_err;
}

static int gt_rdb_stmt_sqlite_bind_string(GtRDBStmt *st, unsigned long param_no,
                                          const char *val, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int rval, had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  gt_assert(param_no < sts->num_params);
  rval = sqlite3_bind_text(sts->stmt, param_no + 1, val, strlen(val),
                           SQLITE_STATIC);
  if (rval != SQLITE_OK) {
    gt_error_set(err, GT_SQLITE_ERRMSG, rval, sqlite3_errmsg(sts->db));
    had_err = -1;
  }
  return had_err;
}

static int gt_rdb_stmt_sqlite_bind_double(GtRDBStmt *st, unsigned long param_no,
                                          double val, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int rval, had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  gt_assert(param_no < sts->num_params);
  rval = sqlite3_bind_double(sts->stmt, param_no + 1, val);
  if (rval != SQLITE_OK) {
    gt_error_set(err, GT_SQLITE_ERRMSG, rval, sqlite3_errmsg(sts->db));
    had_err = -1;
  }
  return had_err;
}

static int gt_rdb_stmt_sqlite_exec(GtRDBStmt *st, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int rval, had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  switch ((rval = sqlite3_step(sts->stmt))) {
    case SQLITE_ROW:
      break;
    case SQLITE_DONE:
      had_err = 1;  /* last row read */
      break;
    default:
      gt_error_set(err, GT_SQLITE_ERRMSG, rval, sqlite3_errmsg(sts->db));
      had_err = -1;
      break;
  }
  return had_err;
}

static int gt_rdb_stmt_sqlite_get_int(GtRDBStmt *st, unsigned long field_no,
                                      int *result, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  if (!sts->stmt) {
    gt_error_set(err, "uninitialized statement");
    had_err = -1;
  }
  if (!had_err) {
    *result = sqlite3_column_int(sts->stmt, field_no);
  }
  return had_err;
}

static int gt_rdb_stmt_sqlite_get_ulong(GtRDBStmt *st, unsigned long field_no,
                                        unsigned long *result, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  if (!sts->stmt) {
    gt_error_set(err, "uninitialized statement");
    had_err = -1;
  }
  if (!had_err) {
    *result = (unsigned long) sqlite3_column_int(sts->stmt, field_no);
  }
  return had_err;
}

static int gt_rdb_stmt_sqlite_get_string(GtRDBStmt *st, unsigned long field_no,
                                         GtStr *result, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  if (!sts->stmt) {
    gt_error_set(err, "uninitialized statement");
    had_err = -1;
  }
  if (!had_err) {
    gt_str_reset(result);
    gt_str_append_cstr(result,
                       (const char*) sqlite3_column_text(sts->stmt, field_no));
  }
  return had_err;
}

static int gt_rdb_stmt_sqlite_get_double(GtRDBStmt *st, unsigned long field_no,
                                         double *result, GtError *err)
{
  GtRDBStmtSqlite *sts;
  int had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  sts = gt_rdb_stmt_sqlite_cast(st);
  if (!sts->stmt) {
    gt_error_set(err, "uninitialized statement");
    had_err = -1;
  }
  if (!had_err) {
    *result = sqlite3_column_double(sts->stmt, field_no);
  }
  return had_err;
}

static GtCstrTable* gt_rdb_sqlite_get_indexes(GtRDB *rdb, GtError *err)
{
  GtRDBStmt *stmt;
  GtCstrTable *tab;
  int rval = 0;
  gt_assert(rdb);
  gt_error_check(err);
  if ((stmt = gt_rdb_prepare(rdb, "SELECT name FROM sqlite_master "
                                  "WHERE type='index'",
                                  0,
                                  err)) == NULL)
    return NULL;
  tab = gt_cstr_table_new();
  while (!rval) {
    GtStr *key;
    rval = gt_rdb_stmt_exec(stmt, err);
    if (rval) break;
    key = gt_str_new();
    gt_rdb_stmt_get_string(stmt, 0, key, err);
    gt_cstr_table_add(tab, gt_str_get(key));
    gt_str_delete(key);
  }
  if (rval < 0) {
    gt_cstr_table_delete(tab);
    gt_rdb_stmt_delete(stmt);
    return NULL;
  }
  gt_rdb_stmt_delete(stmt);
  return tab;
}

static GtCstrTable* gt_rdb_sqlite_get_tables(GtRDB *rdb, GtError *err)
{
  GtRDBStmt *stmt;
  GtCstrTable *tab;
  int rval = 0;
  gt_assert(rdb);
  gt_error_check(err);
  if ((stmt = gt_rdb_prepare(rdb, "SELECT name FROM sqlite_master "
                                  "WHERE type='table'",
                                  0,
                                  err)) == NULL)
    return NULL;
  tab = gt_cstr_table_new();
  while (!rval) {
    GtStr *key;
    rval = gt_rdb_stmt_exec(stmt, err);
    if (rval) break;
    key = gt_str_new();
    gt_rdb_stmt_get_string(stmt, 0, key, err);
    gt_cstr_table_add(tab, gt_str_get(key));
    gt_str_delete(key);
  }
  if (rval < 0) {
    gt_cstr_table_delete(tab);
    gt_rdb_stmt_delete(stmt);
    return NULL;
  }
  gt_rdb_stmt_delete(stmt);
  return tab;
}

const GtRDBClass* gt_rdb_sqlite_class(void)
{
  static const GtRDBClass *rdbs = NULL;
  gt_class_alloc_lock_enter();
  if (!rdbs) {
    rdbs = gt_rdb_class_new(sizeof (GtRDBSqlite),
                            gt_rdb_sqlite_delete,
                            gt_rdb_sqlite_prepare,
                            gt_rdb_sqlite_last_inserted_id,
                            gt_rdb_sqlite_accept,
                            gt_rdb_sqlite_get_indexes,
                            gt_rdb_sqlite_get_tables);
  }
  gt_class_alloc_lock_leave();
  return rdbs;
}

const GtRDBStmtClass* gt_rdb_stmt_sqlite_class(void)
{
  static const GtRDBStmtClass *rdbss = NULL;
  gt_class_alloc_lock_enter();
  if (!rdbss) {
    rdbss = gt_rdb_stmt_class_new(sizeof (GtRDBStmtSqlite),
                                  gt_rdb_stmt_sqlite_reset,
                                  gt_rdb_stmt_sqlite_bind_int,
                                  gt_rdb_stmt_sqlite_bind_ulong,
                                  gt_rdb_stmt_sqlite_bind_string,
                                  gt_rdb_stmt_sqlite_bind_double,
                                  gt_rdb_stmt_sqlite_exec,
                                  gt_rdb_stmt_sqlite_get_int,
                                  gt_rdb_stmt_sqlite_get_ulong,
                                  gt_rdb_stmt_sqlite_get_string,
                                  gt_rdb_stmt_sqlite_get_double,
                                  gt_rdb_stmt_sqlite_delete);
  }
  gt_class_alloc_lock_leave();
  return rdbss;
}

#endif
