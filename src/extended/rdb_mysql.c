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

#ifdef HAVE_MYSQL
#include <string.h>
#include <mysql/mysql.h>
#include "core/class_alloc_lock.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/hashtable.h"
#include "core/str_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/rdb_mysql_api.h"
#include "extended/rdb_rep.h"

struct GtRDBMySQL {
  const GtRDB parent_instance;
  MYSQL conn;
  GtStr *database;
};

struct GtRDBStmtMySQL {
    const GtRDBStmt parent_instance;
    unsigned long num_params;
    GtHashtable *buffers, *returned_strings;
    bool executed;
    my_bool update_maxlengths;
    GtStr *query;
    MYSQL_BIND *params, *results;
    MYSQL_STMT *stmt;
    MYSQL *conn;
};

const GtRDBClass*     gt_rdb_mysql_class(void);
const GtRDBStmtClass* gt_rdb_stmt_mysql_class(void);

#define GT_MYSQL_ERRMSG  "MySQL error code %d: %s"

#define gt_rdb_mysql_cast(RDB)\
        gt_rdb_cast(gt_rdb_mysql_class(), RDB)

#define gt_rdb_stmt_mysql_cast(RDB)\
        gt_rdb_stmt_cast(gt_rdb_stmt_mysql_class(), RDB)

GtRDB* gt_rdb_mysql_new(const char *server,
                        unsigned int port,
                        const char *database,
                        const char *username,
                        const char *password,
                        GtError *err)
{
  GtRDBMySQL *rdbm;
  GtRDB *rdb = gt_rdb_create(gt_rdb_mysql_class());
  rdbm = gt_rdb_mysql_cast(rdb);
  mysql_init(&rdbm->conn);
  mysql_options(&rdbm->conn, MYSQL_READ_DEFAULT_GROUP, "genometools");
  if (!mysql_real_connect(&rdbm->conn, server, username, password, database,
                          port, NULL,
                          CLIENT_COMPRESS | CLIENT_MULTI_STATEMENTS))
  {
    gt_error_set(err, "cannot connect to database: %s",
                 mysql_error(&rdbm->conn));
    mysql_close(&rdbm->conn);
    gt_rdb_delete(rdb);
    return NULL;
  }
  rdbm->database = gt_str_new_cstr(database);
  return rdb;
}

static void free_buf(void *elem)
{
  gt_free(*(void**) elem);
}

static void free_str(void *elem)
{
  gt_str_delete(*(GtStr**) elem);
}

static GtRDBStmt* gt_rdb_mysql_prepare(GtRDB *rdb, const char *query,
                                       unsigned long num_params, GtError *err)
{
  GtRDBStmt *st = NULL;
  GtRDBStmtMySQL *stm = NULL;
  GtRDBMySQL *rdbm;
  int had_err = 0, retval = 0;
  /* we need these to keep track of result/parameter and string buffers */
  HashElemInfo str_buffer_hash = {
      gt_ht_ptr_elem_hash,
      { free_str },
      sizeof (GtStr*),
      gt_ht_ptr_elem_cmp,
      NULL,
      NULL
    },
    buffer_hash = {
      gt_ht_ptr_elem_hash,
      { free_buf },
      sizeof (void*),
      gt_ht_ptr_elem_cmp,
      NULL,
      NULL
    };
  MYSQL_STMT *tmp = NULL;
  gt_assert(rdb && query);
  gt_error_check(err);

  rdbm = gt_rdb_mysql_cast(rdb);
  tmp = mysql_stmt_init(&rdbm->conn);
  if ((retval = mysql_stmt_prepare(tmp, query, strlen(query)))) {
    gt_error_set(err, GT_MYSQL_ERRMSG, retval, mysql_stmt_error(tmp));
    had_err = -1;
  }
  if (!had_err) {
    int param_count;
    param_count = mysql_stmt_param_count(tmp);
    if (param_count != num_params) {
      gt_error_set(err, "invalid parameter count: %lu expected, %d given",
                   num_params, param_count);
      mysql_stmt_close(tmp);
      had_err = -1;
    }
  }
  if (!had_err) {
    st = gt_rdb_stmt_create(gt_rdb_stmt_mysql_class());
    stm = gt_rdb_stmt_mysql_cast(st);
    stm->num_params = num_params;
    stm->query = gt_str_new_cstr(query);
    stm->buffers = gt_hashtable_new(buffer_hash);
    stm->returned_strings = gt_hashtable_new(str_buffer_hash);
    stm->stmt = tmp;
    stm->update_maxlengths = true;
    stm->params = gt_calloc(num_params, sizeof (MYSQL_BIND));
    mysql_stmt_attr_set(tmp, STMT_ATTR_UPDATE_MAX_LENGTH,
                        &stm->update_maxlengths);
    memset(stm->params, 0, num_params*sizeof (MYSQL_BIND));
    stm->conn = &rdbm->conn;
  }
  return st;
}

static unsigned long gt_rdb_mysql_last_inserted_id(GtRDB *rdb,
                                                   GT_UNUSED const char *table,
                                                   GT_UNUSED GtError *err)
{
  GtRDBMySQL *rdbm;
  gt_assert(rdb);
  gt_error_check(err);
  rdbm = gt_rdb_mysql_cast(rdb);
  return mysql_insert_id(&rdbm->conn); /* TODO: find out whether this is better
                                                replaced by INFORMATION_SCHEMA
                                                query */
}

static GtCstrTable* gt_rdb_mysql_get_tables(GtRDB *rdb, GtError *err)
{
  GtRDBMySQL *rdbm;
  MYSQL_RES *res;
  MYSQL_ROW row;
  GtCstrTable *tab;
  gt_assert(rdb);
  gt_error_check(err);
  rdbm = gt_rdb_mysql_cast(rdb);
  gt_assert(&rdbm->conn);
  res = mysql_list_tables(&rdbm->conn, NULL); /* NULL means 'all tables' */
  if (!res) {
    gt_error_set(err, "error trying to list tables: %s",
                 mysql_error(&rdbm->conn));
    return NULL;
  }
  tab = gt_cstr_table_new();
  while ((row = mysql_fetch_row(res))) {
    char buf[BUFSIZ];
    unsigned long *lengths;
    memset(buf, 0, BUFSIZ);
    lengths = mysql_fetch_lengths(res);
    (void) snprintf(buf, MIN(BUFSIZ, lengths[0])*sizeof (char), "%s",
                    (char*) row[0] ? (char*) row[0] : "NULL");
    gt_cstr_table_add(tab, buf);
  }
  mysql_free_result(res);
  return tab;
}

static GtCstrTable* gt_rdb_mysql_get_indexes(GtRDB *rdb, GtError *err)
{
  GT_UNUSED GtRDBMySQL *rdbm;
  GtRDBStmt *stmt;
  GtCstrTable *tab;
  int rval = 0;
  gt_assert(rdb);
  gt_error_check(err);
  rdbm = gt_rdb_mysql_cast(rdb);
  /* TODO: implement a way to do this in MySQL 4 */
  if ((stmt = gt_rdb_prepare(rdb, "SELECT DISTINCT INDEX_NAME "
                                  "FROM INFORMATION_SCHEMA.STATISTICS ",
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

static int gt_rdb_mysql_accept(GtRDB *rdb, GtRDBVisitor *v, GtError *err)
{
  GtRDBMySQL *rdbm;
  gt_assert(rdb && v);
  gt_error_check(err);
  rdbm = gt_rdb_mysql_cast(rdb);
  return gt_rdb_visitor_visit_mysql(v, rdbm, err);
}

static void gt_rdb_mysql_delete(GtRDB *rdb)
{
  GtRDBMySQL *rdbm;
  if (!rdb) return;
  rdbm = gt_rdb_mysql_cast(rdb);
  mysql_close(&rdbm->conn);
  gt_str_delete(rdbm->database);
}

static void gt_rdb_stmt_mysql_delete(GtRDBStmt *st)
{
  GtRDBStmtMySQL *stm;
  if (!st) return;
  stm = gt_rdb_stmt_mysql_cast(st);
  if (stm->stmt) {
    mysql_stmt_free_result(stm->stmt);
    mysql_stmt_close(stm->stmt);
  }
  gt_free(stm->params);
  gt_free(stm->results);
  gt_hashtable_delete(stm->buffers);
  gt_hashtable_delete(stm->returned_strings);
  gt_str_delete(stm->query);
}

static int gt_rdb_stmt_mysql_reset(GtRDBStmt *st, GtError *err)
{
  GtRDBStmtMySQL *stm;
  int rval, had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  gt_hashtable_reset(stm->buffers);
  gt_hashtable_reset(stm->returned_strings);
  mysql_stmt_free_result(stm->stmt);
  if ((rval = mysql_stmt_reset(stm->stmt))) {
    gt_error_set(err, GT_MYSQL_ERRMSG, rval, mysql_stmt_error(stm->stmt));
    had_err = -1;
  }
  memset(stm->params, 0, stm->num_params*sizeof (MYSQL_BIND));
  gt_free(stm->results);
  stm->results = NULL;
  if (!had_err)
    stm->executed = false;
  return had_err;
}

static int gt_rdb_stmt_mysql_bind_int(GtRDBStmt *st, unsigned long param_no,
                                      int val, GT_UNUSED GtError *err)
{
  GtRDBStmtMySQL *stm;
  int *lval, had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  gt_assert(param_no < stm->num_params);
  lval = gt_malloc(sizeof (int));
  *lval = val;
  gt_hashtable_add(stm->buffers, &lval);
  stm->params[param_no].buffer_type = MYSQL_TYPE_LONG;
  stm->params[param_no].buffer = (char*) lval;
  stm->params[param_no].is_null = 0;
  stm->params[param_no].length = 0;
  return had_err;
}

static int gt_rdb_stmt_mysql_bind_ulong(GtRDBStmt *st, unsigned long param_no,
                                        unsigned long val,
                                        GT_UNUSED GtError *err)
{
  GtRDBStmtMySQL *stm;
  unsigned long *lval;
  int had_err = 0;
  gt_assert(st);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  gt_assert(param_no < stm->num_params);
  lval = gt_malloc(sizeof (unsigned long));
  *lval = val;
  gt_hashtable_add(stm->buffers, &lval);
  stm->params[param_no].buffer_type = MYSQL_TYPE_LONG;
  stm->params[param_no].is_unsigned = true;
  stm->params[param_no].buffer = (char*) lval;
  stm->params[param_no].is_null = 0;
  stm->params[param_no].length = 0;
  return had_err;
}

static int gt_rdb_stmt_mysql_bind_string(GtRDBStmt *st, unsigned long param_no,
                                         const char *val,
                                         GT_UNUSED GtError *err)
{
  GtRDBStmtMySQL *stm;
  int had_err = 0;
  char *str;
  unsigned long *length;
  gt_assert(st);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  gt_assert(param_no < stm->num_params);
  /* allocate buffer for string until execution */
  str = gt_calloc(strlen(val)+1, sizeof (char));
  strncpy(str, val, strlen(val));
  gt_hashtable_add(stm->buffers, &str);
  /* allocate space for length */
  length = gt_malloc(sizeof (unsigned long));
  *length = strlen(str);
  gt_hashtable_add(stm->buffers, &length);
  /* fill param structure */
  stm->params[param_no].buffer_type = MYSQL_TYPE_STRING;
  stm->params[param_no].buffer = str;
  stm->params[param_no].is_null = 0;
  stm->params[param_no].buffer_length = strlen(str);
  stm->params[param_no].length = length;
  return had_err;
}

static int gt_rdb_stmt_mysql_bind_double(GtRDBStmt *st, unsigned long param_no,
                                         double val, GT_UNUSED GtError *err)
{
  GtRDBStmtMySQL *stm;
  int had_err = 0;
  double *lval;
  gt_assert(st);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  gt_assert(param_no < stm->num_params);
  lval = gt_malloc(sizeof (double));
  *lval = val;
  gt_hashtable_add(stm->buffers, &lval);
  stm->params[param_no].buffer_type = MYSQL_TYPE_DOUBLE;
  stm->params[param_no].is_unsigned = false;
  stm->params[param_no].buffer = (char*) lval;
  stm->params[param_no].is_null = 0;
  stm->params[param_no].length = 0;
  return had_err;
}

#define CHECK_INIT_STATEMENT \
  if (!stm->stmt || !stm->executed) { \
    gt_error_set(err, "uninitialized statement"); \
    had_err = -1; \
  } \
  if (!had_err && !stm->results) { \
    gt_error_set(err, "results not fetched!"); \
    had_err = -1; \
  }

static int gt_rdb_stmt_mysql_get_int(GtRDBStmt *st, unsigned long field_no,
                                     int *result, GtError *err)
{
  GtRDBStmtMySQL *stm;
  int had_err = 0;
  gt_assert(st && result);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  CHECK_INIT_STATEMENT
  if (!had_err       /* TODO: check these VVV */
       && stm->results[field_no].buffer_type != MYSQL_TYPE_LONG
       && stm->results[field_no].buffer_type != MYSQL_TYPE_TINY
       && stm->results[field_no].buffer_type != MYSQL_TYPE_INT24
       && stm->results[field_no].buffer_type != MYSQL_TYPE_SHORT) {
    gt_error_set(err, "incompatible type: %d!",
                 stm->results[field_no].buffer_type);
    had_err = -1;
  }
  if (!had_err) {
    switch (stm->results[field_no].buffer_type) {
      case MYSQL_TYPE_LONG:
      case MYSQL_TYPE_INT24:
        {int val = *(int*) stm->results[field_no].buffer;
        *result = val;}
        break;
      case MYSQL_TYPE_SHORT:
        {short int val = *(short int*) stm->results[field_no].buffer;
        *result = (int) val;}
        break;
      case MYSQL_TYPE_TINY:
        {signed char val = *(signed char*) stm->results[field_no].buffer;
        *result = (int) val;}
        break;
      default:
        gt_assert(false); /* should not happen */
    }
  }
  return had_err;
}

static int gt_rdb_stmt_mysql_get_ulong(GtRDBStmt *st, unsigned long field_no,
                                       unsigned long *result, GtError *err)
{
  GtRDBStmtMySQL *stm;
  int had_err = 0;
  gt_assert(st && result);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  CHECK_INIT_STATEMENT
  if (!had_err && stm->results[field_no].buffer_type != MYSQL_TYPE_LONG) {
    gt_error_set(err, "incompatible type!");
    had_err = -1;
  }
  if (!had_err)
    *result = *(unsigned long*) stm->results[field_no].buffer;
  return had_err;
}

static int gt_rdb_stmt_mysql_get_string(GtRDBStmt *st, unsigned long field_no,
                                        GtStr *result, GtError *err)
{
  GtRDBStmtMySQL *stm;
  int had_err = 0;
  gt_assert(st && result);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  CHECK_INIT_STATEMENT
  if (!had_err
        && stm->results[field_no].buffer_type != MYSQL_TYPE_STRING
        && stm->results[field_no].buffer_type != MYSQL_TYPE_VAR_STRING
        && stm->results[field_no].buffer_type != MYSQL_TYPE_BLOB
        && stm->results[field_no].buffer_type != MYSQL_TYPE_TINY_BLOB
        && stm->results[field_no].buffer_type != MYSQL_TYPE_MEDIUM_BLOB
        && stm->results[field_no].buffer_type != MYSQL_TYPE_LONG_BLOB
        && stm->results[field_no].buffer_type != MYSQL_TYPE_BIT)
  {
    gt_error_set(err, "incompatible type!");
    had_err = -1;
  }
  if (!had_err) {
    gt_str_reset(result);
    gt_str_append_cstr_nt(result,
                          (char*)stm->results[field_no].buffer,
                          *stm->results[field_no].length);
  }
  return had_err;
}

static int gt_rdb_stmt_mysql_get_double(GtRDBStmt *st, unsigned long field_no,
                                        double *result, GtError *err)
{
  GtRDBStmtMySQL *stm;
  MYSQL_BIND res_bind;
  int had_err = 0;
  double res_double = GT_UNDEF_DOUBLE;
  my_bool error, is_null;
  gt_assert(st && result);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  CHECK_INIT_STATEMENT
  if (!had_err) {
    memset(&res_bind, 0, sizeof (res_bind));
    res_bind.buffer_type = MYSQL_TYPE_DOUBLE;
    res_bind.buffer = &res_double;
    res_bind.error = &error;
    res_bind.is_null = &is_null;
    if ((had_err = mysql_stmt_fetch_column(stm->stmt, &res_bind, field_no, 0)))
      gt_error_set(err, GT_MYSQL_ERRMSG, had_err, mysql_stmt_error(stm->stmt));
  }
  if (!had_err)
    *result = res_double;
  return had_err;
}

static int gt_rdb_stmt_mysql_exec(GtRDBStmt *st, GtError *err)
{
  GtRDBStmtMySQL *stm;
  int rval, had_err = 0, num_fields;
  MYSQL_RES *meta_res = NULL;
  gt_assert(st);
  gt_error_check(err);
  stm = gt_rdb_stmt_mysql_cast(st);
  if (!stm->executed) {
    if (stm->num_params > 0) {
      gt_assert(stm->stmt && stm->params);
      if ((rval = mysql_stmt_bind_param(stm->stmt, stm->params))) {
        gt_error_set(err, GT_MYSQL_ERRMSG, rval, mysql_stmt_error(stm->stmt));
        had_err = -1;
      }
    }
    if (!had_err && (rval = mysql_stmt_execute(stm->stmt))) {
      gt_error_set(err, GT_MYSQL_ERRMSG, rval, mysql_stmt_error(stm->stmt));
      had_err = -1;
    }
    if (!had_err) {
      stm->executed = true;
      if (mysql_stmt_store_result(stm->stmt)) {
        gt_error_set(err, GT_MYSQL_ERRMSG,
                     had_err, mysql_stmt_error(stm->stmt));
        had_err = -1;
      }
      meta_res = mysql_stmt_result_metadata(stm->stmt);
      if (!had_err && meta_res) {
        int i = 0;
        /* statement returned a result */
        num_fields = mysql_num_fields(meta_res);
        stm->results = gt_calloc(num_fields, sizeof (MYSQL_BIND));
        /* prepare result buffers for each field */
        for (i=0;i<num_fields;i++) {
          MYSQL_FIELD *field;
          field = mysql_fetch_field(meta_res);
          stm->results[i].buffer_type = field->type;
          switch (field->type) {
            case MYSQL_TYPE_DOUBLE:
              {double *dbl = gt_calloc(1, sizeof (double));
              gt_hashtable_add(stm->buffers, &dbl);
              stm->results[i].buffer_length = sizeof (double);
              stm->results[i].buffer = dbl;}
              break;
            case MYSQL_TYPE_LONG:
            case MYSQL_TYPE_INT24:
            {int *l = gt_calloc(1, sizeof (int));
              gt_hashtable_add(stm->buffers, &l);
              stm->results[i].is_unsigned = false;
              stm->results[i].buffer_length = sizeof (int);
              stm->results[i].buffer = l;}
            case MYSQL_TYPE_SHORT:
            {short int *l = gt_calloc(1, sizeof (short int));
              gt_hashtable_add(stm->buffers, &l);
              stm->results[i].is_unsigned = false;
              stm->results[i].buffer_length = sizeof (short int);
              stm->results[i].buffer = l;}
            case MYSQL_TYPE_TINY:
              {signed char *l = gt_calloc(1, sizeof (signed char));
              gt_hashtable_add(stm->buffers, &l);
              stm->results[i].is_unsigned = false;
              stm->results[i].buffer_length = sizeof (signed char);
              stm->results[i].buffer = l;}
              break;
            case MYSQL_TYPE_STRING:
            case MYSQL_TYPE_VAR_STRING:
            case MYSQL_TYPE_BLOB:
            case MYSQL_TYPE_TINY_BLOB:
            case MYSQL_TYPE_MEDIUM_BLOB:
            case MYSQL_TYPE_LONG_BLOB:
            case MYSQL_TYPE_BIT:
              {char *str = gt_calloc(field->max_length+1, sizeof (char));
              gt_hashtable_add(stm->buffers, &str);
              unsigned long *length = gt_calloc(1, sizeof (unsigned long));
              gt_hashtable_add(stm->buffers, &length);
              stm->results[i].buffer = str;
              stm->results[i].buffer_length = field->max_length;
              stm->results[i].length = length;}
              break;
            default:
              /* unsupported data type */
              break;
          }
        }
        if (!had_err)
          mysql_stmt_bind_result(stm->stmt, stm->results);
        mysql_free_result(meta_res);
      } else {
        return 1;
      }
    }
  }
  if (!had_err) {
    switch ((rval = mysql_stmt_fetch(stm->stmt))) {
      case 0:
      default:
        break;
      case MYSQL_NO_DATA:
        had_err = 1;  /* last row read */
        break;
      case 1:
        gt_error_set(err, GT_MYSQL_ERRMSG, mysql_stmt_errno(stm->stmt),
                     mysql_stmt_error(stm->stmt));
        had_err = -1;
        break;
    }
  }
  return had_err;
}

const GtRDBClass* gt_rdb_mysql_class(void)
{
  static const GtRDBClass *rdbm = NULL;
  gt_class_alloc_lock_enter();
  if (!rdbm) {
    rdbm = gt_rdb_class_new(sizeof (GtRDBMySQL),
                            gt_rdb_mysql_delete,
                            gt_rdb_mysql_prepare,
                            gt_rdb_mysql_last_inserted_id,
                            gt_rdb_mysql_accept,
                            gt_rdb_mysql_get_indexes,
                            gt_rdb_mysql_get_tables);
  }
  gt_class_alloc_lock_leave();
  return rdbm;
}

const GtRDBStmtClass* gt_rdb_stmt_mysql_class(void)
{
  static const GtRDBStmtClass *rdbms = NULL;
  gt_class_alloc_lock_enter();
  if (!rdbms) {
    rdbms = gt_rdb_stmt_class_new(sizeof (GtRDBStmtMySQL),
                                  gt_rdb_stmt_mysql_reset,
                                  gt_rdb_stmt_mysql_bind_int,
                                  gt_rdb_stmt_mysql_bind_ulong,
                                  gt_rdb_stmt_mysql_bind_string,
                                  gt_rdb_stmt_mysql_bind_double,
                                  gt_rdb_stmt_mysql_exec,
                                  gt_rdb_stmt_mysql_get_int,
                                  gt_rdb_stmt_mysql_get_ulong,
                                  gt_rdb_stmt_mysql_get_string,
                                  gt_rdb_stmt_mysql_get_double,
                                  gt_rdb_stmt_mysql_delete);
  }
  gt_class_alloc_lock_leave();
  return rdbms;
}

#endif
