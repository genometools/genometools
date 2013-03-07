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

#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/fileutils_api.h"
#include "core/hashtable.h"
#include "core/hashmap.h"
#include "core/hashmap-generic.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/range.h"
#include "core/strand_api.h"
#include "core/thread_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/fa.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "core/xposix.h"
#include "extended/anno_db_gfflike_api.h"
#include "extended/anno_db_prepared_stmt.h"
#include "extended/anno_db_schema_rep.h"
#include "extended/feature_index_rep.h"
#include "extended/feature_index.h"
#include "extended/feature_visitor.h"
#include "extended/feature_node.h"
#include "extended/feature_node_observer.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/rdb_api.h"
#include "extended/rdb_sqlite_api.h"
#include "extended/rdb_visitor_rep.h"

struct GtAnnoDBGFFlike {
  const GtAnnoDBSchema parent_instance;
  GtRDB *db;
  GtRDBVisitor *visitor;
};

typedef struct {
  const GtRDBVisitor parent_instance;
  GtAnnoDBGFFlike *annodb;
} GFFlikeSetupVisitor;

typedef struct {
  const GtFeatureIndex parent_instance;
  GtHashmap *node_to_parent_array,
            *seqid_cache,
            *source_cache,
            *string_caches,
            *ref_nodes,
            *added,
            *deleted,
            *changed;
  GtHashtable *cache_node2id,
              *cache_id2node;
  GtRDBStmt *stmts[GT_PSTMT_NOF_STATEMENTS];
  GtFeatureNodeObserver *obs;
  GtRDB *db;
  GtMutex *dblock;
  bool transaction_lock;
} GtFeatureIndexGFFlike;

const GtAnnoDBSchemaClass* gt_anno_db_gfflike_class(void);
static const GtRDBVisitorClass* gfflike_setup_visitor_class(void);
static const GtFeatureIndexClass* feature_index_gfflike_class(void);

#define anno_db_gfflike_cast(V)\
        gt_anno_db_schema_cast(gt_anno_db_gfflike_class(), V)

#define gfflike_setup_visitor_cast(V)\
        gt_rdb_visitor_cast(gfflike_setup_visitor_class(), V)

#define feature_index_gfflike_cast(V)\
        gt_feature_index_cast(feature_index_gfflike_class(), V)

static int anno_db_gfflike_validate_sqlite(GtRDBSqlite *db, GtError *err,
                                           bool *check)
{
  int had_err = 0;
  bool ok = true;
  GtCstrTable *tables = gt_rdb_get_tables((GtRDB*) db, err);
  if (!tables)
    had_err = -1;
  if (ok && !had_err && !gt_cstr_table_get(tables, "features")) ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "types")) ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "sequenceregions"))
    ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "parents")) ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "sources")) ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "attributes")) ok = false;
  /* TODO: more elaborate checks? */
  if (!had_err)
    *check = ok;
  gt_cstr_table_delete(tables);
  return had_err;
}

static int anno_db_gfflike_create_tables_sqlite(GtRDBSqlite *db, GtError *err)
{
  int had_err = 0;
  GtRDBStmt *stmt;

  stmt = gt_rdb_prepare((GtRDB*) db, "PRAGMA synchronous=OFF", 0, err);
  if (!stmt) return -1;
  gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db, "PRAGMA default_cache_size=256000", 0,
                        err);
  if (!stmt) return -1;
  gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db, "PRAGMA cache_size=512000", 0, err);
  if (!stmt) return -1;
  gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db, "PRAGMA count_changes=OFF", 0, err);
  if (!stmt) return -1;
  gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db, "PRAGMA temp_store=MEMORY", 0, err);
  if (!stmt) return -1;
  gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db, "PRAGMA journal_mode=MEMORY", 0, err);
  if (!stmt) return -1;
  gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS features "
                           "(id INTEGER PRIMARY KEY AUTOINCREMENT, "
                             "seqid INTEGER NOT NULL "
                             "REFERENCES sequenceregions (sequenceregion_id), "
                           "source INTEGER NOT NULL "
                             "REFERENCES sources (source_id), "
                           "type INTEGER NOT NULL "
                             "REFERENCES types (type_id), "
                           "start INTEGER NOT NULL, "
                           "end INTEGER NOT NULL, "
                           "score REAL NOT NULL, "
                           "strand VARCHAR(1) NOT NULL, "
                           "phase INTEGER NOT NULL, "
                           "is_multi INTEGER NOT NULL, "
                           "is_pseudo INTEGER NOT NULL, "
                           "is_marked INTEGER NOT NULL, "
                           "multi_representative INTEGER NOT NULL)",
                           0, err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS types "
                           "(type_id INTEGER PRIMARY KEY AUTOINCREMENT, "
                           "type_name VARCHAR(255))",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS parents "
                           "(feature_id INTEGER "
                             "REFERENCES features (id), "
                           "parent INTEGER "
                             "REFERENCES features (id), "
                           "PRIMARY KEY (feature_id, parent))",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS sources "
                           "(source_id INTEGER PRIMARY KEY AUTOINCREMENT, "
                           "source_name VARCHAR(255))",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS sequenceregions "
                           "(sequenceregion_id INTEGER "
                             "PRIMARY KEY AUTOINCREMENT, "
                           "sequenceregion_name VARCHAR(255), "
                           "start INTEGER, "
                           "stop INTEGER)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS attributes "
                           "(feature_id INTEGER "
                             "REFERENCES features (id), "
                           "keystr VARCHAR(255), "
                           "value VARCHAR(255), "
                           "PRIMARY KEY(feature_id, keystr))",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  return 0;
}

static int anno_db_gfflike_create_indexes_sqlite(GtRDBSqlite *db, GtError *err)
{
  int had_err = 0;
  GtRDBStmt *stmt;
  gt_assert(db);

  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS feature_all "
                           "ON features (id, start, end, seqid, source, type)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS name_type "
                           "ON types (type_name)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS name_source "
                           "ON sources (source_name)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS feature_seqid "
                           "ON features (seqid)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS name_sequenceregion "
                           "ON sequenceregions (sequenceregion_name)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS attribs_value "
                           "ON attributes (value)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS attribs_key "
                           "ON attributes (keystr)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS attribs_feature "
                           "ON attributes (feature_id)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX IF NOT EXISTS parent_id "
                           "ON parents (feature_id)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  return 0;
}

static int anno_db_gfflike_validate_mysql(GtRDBMySQL *db, GtError *err,
                                          bool *check)
{
  int had_err = 0;
  bool ok = true;
  GtCstrTable *tables = gt_rdb_get_tables((GtRDB*) db, err);
  if (!tables)
    had_err = -1;
  /* XXX: for some reason, MySQL omits the last letter of the table name?! */
  if (ok && !had_err && !gt_cstr_table_get(tables, "feature")) ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "type")) ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "sequenceregion"))
    ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "parent")) ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "source")) ok = false;
  if (ok && !had_err && !gt_cstr_table_get(tables, "attribute")) ok = false;
  /* TODO: more elaborate checks? */
  if (!had_err)
    *check = ok;
  gt_cstr_table_delete(tables);
  return had_err;
}

static int anno_db_gfflike_create_tables_mysql(GtRDBMySQL *db, GtError *err)
{
  int had_err = 0;
  GtRDBStmt *stmt;
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS features "
                           "(id INTEGER AUTO_INCREMENT PRIMARY KEY, "
                             "seqid INTEGER NOT NULL "
                             "REFERENCES sequenceregions (sequenceregion_id), "
                           "source INTEGER NOT NULL "
                             "REFERENCES sources (source_id), "
                           "type INTEGER NOT NULL "
                             "REFERENCES types (type_id), "
                           "start INTEGER NOT NULL, "
                           "end INTEGER NOT NULL, "
                           "score REAL NOT NULL, "
                           "strand INTEGER DEFAULT 3 NOT NULL, "
                           "phase INTEGER NOT NULL, "
                           "is_multi INTEGER NOT NULL, "
                           "is_pseudo INTEGER NOT NULL, "
                           "is_marked INTEGER NOT NULL, "
                           "multi_representative INTEGER NOT NULL)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS types "
                           "(type_id INTEGER AUTO_INCREMENT PRIMARY KEY, "
                           "type_name VARCHAR(255))",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS parents "
                           "(feature_id INTEGER "
                             "REFERENCES features (id), "
                           "parent INTEGER "
                             "REFERENCES features (id), "
                           "PRIMARY KEY (feature_id, parent))",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS sources "
                           "(source_id INTEGER AUTO_INCREMENT PRIMARY KEY, "
                           "source_name VARCHAR(255))",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS sequenceregions "
                           "(sequenceregion_id INTEGER "
                             "AUTO_INCREMENT PRIMARY KEY, "
                           "sequenceregion_name VARCHAR(255), "
                           "start INTEGER, "
                           "stop INTEGER)",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE TABLE IF NOT EXISTS attributes "
                           "(id INTEGER AUTO_INCREMENT PRIMARY KEY, "
                             "feature_id INTEGER "
                             "REFERENCES features (id), "
                           "keystr VARCHAR(255), "
                           "value VARCHAR(255))",
                           0,
                           err);
  if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
    return -1;
  } else gt_rdb_stmt_delete(stmt);
  return 0;
}

static int anno_db_gfflike_create_indexes_mysql(GtRDBMySQL *db, GtError *err)
{
  int had_err = 0;
  GtCstrTable *cst;
  GtRDBStmt *stmt;
  gt_assert(db);

  if (!(cst = gt_rdb_get_indexes((GtRDB*) db, err))) {
    return -1;
  }

  if (!gt_cstr_table_get(cst, "feature_all")) {
    stmt = gt_rdb_prepare((GtRDB*) db,
                             "CREATE INDEX feature_all "
                             "ON features (start, end, seqid, source, type)",
                             0,
                             err);
    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }

  if (!gt_cstr_table_get(cst, "name_type")) {
    stmt = gt_rdb_prepare((GtRDB*) db,
                             "CREATE INDEX name_type "
                             "ON types (type_name)",
                             0,
                             err);
    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }

  if (!gt_cstr_table_get(cst, "name_source")) {
    stmt = gt_rdb_prepare((GtRDB*) db,
                             "CREATE INDEX name_source "
                             "ON sources (source_name)",
                             0,
                             err);
    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }

  if (!gt_cstr_table_get(cst, "feature_seqid")) {
    stmt = gt_rdb_prepare((GtRDB*) db,
                             "CREATE INDEX feature_seqid "
                             "ON features (seqid)",
                             0,
                             err);
    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }

  if (!gt_cstr_table_get(cst, "name_sequenceregion")) {
    stmt = gt_rdb_prepare((GtRDB*) db,
                             "CREATE INDEX name_sequenceregion "
                             "ON sequenceregions (sequenceregion_name)",
                             0,
                             err);
    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }

  if (!gt_cstr_table_get(cst, "attribs_value")) {
    stmt = gt_rdb_prepare((GtRDB*) db,
                             "CREATE INDEX attribs_value "
                             "ON attributes (value)",
                             0,
                             err);

    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }

  if (!gt_cstr_table_get(cst, "attribs_key")) {
    stmt = gt_rdb_prepare((GtRDB*) db,
                             "CREATE INDEX attribs_key "
                             "ON attributes (keystr)",
                             0,
                             err);
    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }

  if (!gt_cstr_table_get(cst, "attribs_feature")) {
    stmt = gt_rdb_prepare((GtRDB*) db,
                             "CREATE INDEX attribs_feature "
                             "ON attributes (feature_id)",
                             0,
                             err);
    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }

  if (!gt_cstr_table_get(cst, "parent_id")) {
  stmt = gt_rdb_prepare((GtRDB*) db,
                           "CREATE INDEX parent_id "
                           "ON parents (feature_id)",
                           0,
                           err);
    if (!stmt || (had_err = gt_rdb_stmt_exec(stmt, err)) < 0) {
      gt_rdb_stmt_delete(stmt);
      gt_cstr_table_delete(cst);
      return -1;
    }
    gt_rdb_stmt_delete(stmt);
  }
  return 0;
}

int anno_db_gfflike_init_sqlite(GT_UNUSED GtRDBVisitor *rdbv, GtRDBSqlite *db,
                                GtError *err)
{
  GtCstrTable *cst = NULL;
  GtStrArray *arr = NULL;
  bool check = true;
  int had_err = 0;
  gt_assert(db);

  cst = gt_rdb_get_tables((GtRDB*) db, err);
  if (!cst) {
    had_err = -1;
  }
  if (!had_err) {
    arr = gt_cstr_table_get_all(cst);
    if (!arr)
      had_err = -1;
  }
  if (!had_err) {
    if (gt_str_array_size(arr) == 0) {
      had_err = anno_db_gfflike_create_tables_sqlite(db, err);
    }
  }
  gt_cstr_table_delete(cst);
  gt_str_array_delete(arr);
  if (!had_err)
    had_err = anno_db_gfflike_validate_sqlite(db, err, &check);
  if (!had_err && !check) {
    gt_error_set(err, "possible corruption in database file: "
                      "tables are missing");
    had_err = -1;
  }
  if (!had_err) {
    had_err = anno_db_gfflike_create_indexes_sqlite(db, err);
  }

  return had_err;
}

int anno_db_gfflike_init_mysql(GT_UNUSED GtRDBVisitor *rdbv, GtRDBMySQL *db,
                               GtError *err)
{
  GtCstrTable *cst = NULL;
  GtStrArray *arr = NULL;
  bool check = true;
  int had_err = 0;
  gt_assert(db);

  cst = gt_rdb_get_tables((GtRDB*) db, err);
  if (!cst) {
    had_err = -1;
  }
  if (!had_err) {
    arr = gt_cstr_table_get_all(cst);
    if (!arr)
      had_err = -1;
  }
  if (!had_err) {
    if (gt_str_array_size(arr) == 0) {
      had_err = anno_db_gfflike_create_tables_mysql(db, err);
    }
  }
  gt_cstr_table_delete(cst);
  gt_str_array_delete(arr);
  if (!had_err)
    had_err = anno_db_gfflike_validate_mysql(db, err, &check);
  if (!had_err && !check) {
    gt_error_set(err, "corrupt database schema: tables are missing");
    had_err = -1;
  }
  if (!had_err) {
    had_err = anno_db_gfflike_create_indexes_mysql(db, err);
  }

  return had_err;
}

void anno_db_gfflike_free(GtAnnoDBSchema *s)
{
  GtAnnoDBGFFlike *adg = anno_db_gfflike_cast(s);
  gt_rdb_visitor_delete(adg->visitor);
}

/* for the node->id cache */
DECLARE_HASHMAP(GtFeatureNode*, node, unsigned long, ul, static, inline)
DEFINE_HASHMAP(GtFeatureNode*, node, unsigned long, ul, gt_ht_ptr_elem_hash,
               gt_ht_ptr_elem_cmp, NULL_DESTRUCTOR, NULL_DESTRUCTOR, static,
               inline)

/* for the id->node cache */
DECLARE_HASHMAP(unsigned long, ul, GtFeatureNode*, node, static, inline)
DEFINE_HASHMAP(unsigned long, ul, GtFeatureNode*, node, gt_ht_ul_elem_hash,
               gt_ht_ul_elem_cmp, NULL_DESTRUCTOR, NULL_DESTRUCTOR, static,
               inline)

int gt_feature_index_gfflike_add_region_node(GtFeatureIndex *gfi,
                                             GtRegionNode *rn,
                                             GtError *err)
{
  char *seqid;
  GtFeatureIndexGFFlike *fi;
  GtRange rng;
  int had_err = 0;
  fi = feature_index_gfflike_cast(gfi);
  gt_assert(fi && rn);
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) rn));
  rng = gt_genome_node_get_range((GtGenomeNode*) rn);
  gt_rdb_stmt_reset(fi->stmts[GT_PSTMT_SEQUENCEREGION_INSERT], err);
  gt_rdb_stmt_bind_string(fi->stmts[GT_PSTMT_SEQUENCEREGION_INSERT],
                          0, seqid, err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_SEQUENCEREGION_INSERT],
                       1, (int) rng.start, err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_SEQUENCEREGION_INSERT],
                       2, (int) rng.end, err);
  had_err = (gt_rdb_stmt_exec(fi->stmts[GT_PSTMT_SEQUENCEREGION_INSERT], err)
                              >= 0 ? 0 : -1);
  return had_err;
}

static inline int
gt_feature_index_gfflike_insert_helper(GtFeatureIndexGFFlike *fis,
                                       GtRDBStmt *prepstmt_s,
                                       GtRDBStmt *prepstmt_i,
                                       const char *value,
                                       const char *tabname,
                                       GtError *err)
{
  GtHashmap *cachemap = NULL;
  int *id, rval;
  gt_assert(fis && prepstmt_i && value);

  if (!(cachemap = gt_hashmap_get(fis->string_caches, prepstmt_i))) {
    cachemap = gt_hashmap_new(GT_HASH_STRING, NULL, (GtFree) gt_free_func);
    gt_hashmap_add(fis->string_caches, prepstmt_i, cachemap);
  }
  if (!(id = gt_hashmap_get(cachemap, value))) {
    id = gt_malloc(sizeof (int));
    *id = GT_UNDEF_INT;
    gt_rdb_stmt_reset(prepstmt_s, err);
    gt_rdb_stmt_bind_string(prepstmt_s, 0, value, err);
    rval = gt_rdb_stmt_exec(prepstmt_s, err);
    switch (rval) {
      case 0:
        gt_rdb_stmt_get_int(prepstmt_s, 0, id, err);
        break;
      case 1:
        gt_rdb_stmt_reset(prepstmt_i, err);
        gt_rdb_stmt_bind_string(prepstmt_i, 0, value, err);
        rval = gt_rdb_stmt_exec(prepstmt_i, err);
        if (rval < 0)
          break;
        if (rval == 1)
          *id = (int) gt_rdb_last_inserted_id(fis->db, tabname, err);
        break;
      default:
        gt_error_set(err, "problem executing prepared statement: %d", rval);
        break;
    }
    if (id != NULL && *id != GT_UNDEF_INT) {
      gt_hashmap_add(cachemap, (char*) value, id);
    }
  }
  gt_assert(*id != GT_UNDEF_INT);
  return *id;
}

static int get_parents(GtFeatureNode *gn, void *data, GT_UNUSED GtError *err)
{
  GtHashmap *parentindex = (GtHashmap*) data;
  GtArray *parent_features = NULL;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *fn = NULL;
  gt_assert(gn && parentindex);
  gt_error_check(err);

  fni = gt_feature_node_iterator_new_direct((GtFeatureNode*) gn);
  while ((fn = gt_feature_node_iterator_next(fni))) {
    parent_features = gt_hashmap_get(parentindex, fn);
    if (!parent_features) {
      parent_features = gt_array_new(sizeof (GtFeatureNode*));
      gt_hashmap_add(parentindex, fn, parent_features);
    }
    gt_array_add(parent_features, gn);
  }
  gt_feature_node_iterator_delete(fni);
  return 0;
}

static int set_parents(void *key, void *value, void *data,
                       GtError *err)
{
  int had_err = 0;
  GtFeatureNode *fn = (GtFeatureNode*) key,
                *parent = NULL;
  GtArray *parent_array = (GtArray*) value;
  GtFeatureIndexGFFlike *fi = (GtFeatureIndexGFFlike*) data;
  unsigned long *num, *parent_id;

  num = node_ul_gt_hashmap_get(fi->cache_node2id, fn);
  gt_assert(num);
  if (parent_array && gt_array_size(parent_array)) {
    unsigned long i;
    int rval = 0;
    for (i=0;i<gt_array_size(parent_array);i++) {
      if (rval < 0) {
        had_err = -1;
        break;
      }
      parent = *(GtFeatureNode**) gt_array_get(parent_array, i);
      parent_id = node_ul_gt_hashmap_get(fi->cache_node2id, parent);
      gt_assert(parent_id);
      /* insert parents */
      gt_rdb_stmt_reset(fi->stmts[GT_PSTMT_PARENT_INSERT], err);
      gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_PARENT_INSERT], 0, *num, err);
      gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_PARENT_INSERT], 1, *parent_id,
                           err);
      rval = gt_rdb_stmt_exec(fi->stmts[GT_PSTMT_PARENT_INSERT], err);
    }
  }
  return had_err;
}

static int insert_single_node(GtFeatureIndexGFFlike *fi,
                              unsigned long *id,
                              GtFeatureNode *fn,
                              GtError *err)
{
  GtRange rng;
  GtStrArray *attribs;
  unsigned long i, *myid;
  int type_id           = GT_UNDEF_INT,
      source_id         = GT_UNDEF_INT,
      sequenceregion_id = GT_UNDEF_INT,
      rval = 0,
      had_err = 0,
      multi = 0,
      multi_rep = 0;

  /* do we have the id in the cache already */
  myid = node_ul_gt_hashmap_get(fi->cache_node2id, fn);
  if (myid) {
    *id = *myid;
    return had_err;
  }

  gt_mutex_lock(fi->dblock);

  /* insert details */
  if (!gt_feature_node_is_pseudo(fn)) {
    /* pseudo-features do not have a type */
    type_id =
      gt_feature_index_gfflike_insert_helper(fi,
                                         fi->stmts[GT_PSTMT_TYPE_SELECT],
                                         fi->stmts[GT_PSTMT_TYPE_INSERT],
                                         gt_feature_node_get_type(fn),
                                         "types",
                                         NULL);
  }

  source_id =
    gt_feature_index_gfflike_insert_helper(fi,
                                           fi->stmts[GT_PSTMT_SOURCE_SELECT],
                                           fi->stmts[GT_PSTMT_SOURCE_INSERT],
                                           gt_feature_node_get_source(fn),
                                           "sources",
                                           NULL);

  sequenceregion_id =
    gt_feature_index_gfflike_insert_helper(fi,
                    fi->stmts[GT_PSTMT_SEQUENCEREGION_SELECT],
                    fi->stmts[GT_PSTMT_SEQUENCEREGION_INSERT],
                    gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) fn)),
                    "sequenceregions",
                    NULL);

  /* insert feature */
  rng = gt_genome_node_get_range((GtGenomeNode*) fn);
  gt_rdb_stmt_reset(fi->stmts[GT_PSTMT_FEATURE_INSERT], err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 0,
                       sequenceregion_id, err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 1, source_id, err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 2, type_id, err);
  gt_rdb_stmt_bind_ulong(fi->stmts[GT_PSTMT_FEATURE_INSERT], 3, rng.start,
                         err);
  gt_rdb_stmt_bind_ulong(fi->stmts[GT_PSTMT_FEATURE_INSERT], 4, rng.end, err);
  gt_rdb_stmt_bind_double(fi->stmts[GT_PSTMT_FEATURE_INSERT], 5,
                          gt_feature_node_score_is_defined(fn) ?
                            gt_feature_node_get_score(fn) :
                            GT_UNDEF_DOUBLE,
                          err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 6,
                       gt_feature_node_get_strand(fn), err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 7,
                       gt_feature_node_get_phase(fn), err);

  /* store multi-feature information */
  if (gt_feature_node_is_multi(fn)) {
    GtFeatureNode *rep;
    multi = 1;
    if ((rep = gt_feature_node_get_multi_representative(fn)) != fn) {
      unsigned long *multi_rep_id;
      multi_rep_id = node_ul_gt_hashmap_get(fi->cache_node2id, rep);
      gt_assert(multi_rep_id);
      multi_rep = *multi_rep_id;
    }
  } else multi = 0;

  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 8, multi, err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 9, multi_rep, err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 10,
                       gt_feature_node_is_pseudo(fn), err);
  gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_FEATURE_INSERT], 11,
                       gt_feature_node_is_marked(fn), err);
  rval = gt_rdb_stmt_exec(fi->stmts[GT_PSTMT_FEATURE_INSERT], err);
  if (rval < 0) gt_error_check(err);

  *id = gt_rdb_last_inserted_id(fi->db, "features", err);
  /* cache DB keys to avoid redundant saving of nodes with
     multiple parents */
  node_ul_gt_hashmap_add(fi->cache_node2id, fn, *id);
  ul_node_gt_hashmap_add(fi->cache_id2node, *id, fn);

  /* insert attributes */
  attribs = gt_feature_node_get_attribute_list(fn);
  for (i=0;i<gt_str_array_size(attribs);i++) {
    const char *attr;
    gt_rdb_stmt_reset(fi->stmts[GT_PSTMT_ATTRIBUTE_INSERT], err);
    attr = gt_str_array_get(attribs, i);
    gt_rdb_stmt_bind_int(fi->stmts[GT_PSTMT_ATTRIBUTE_INSERT], 0, *id, err);
    gt_rdb_stmt_bind_string(fi->stmts[GT_PSTMT_ATTRIBUTE_INSERT], 1, attr,
                            err);
    gt_rdb_stmt_bind_string(fi->stmts[GT_PSTMT_ATTRIBUTE_INSERT], 2,
                            gt_feature_node_get_attribute(fn, attr), err);
    rval = gt_rdb_stmt_exec(fi->stmts[GT_PSTMT_ATTRIBUTE_INSERT], err);
    if (rval < 0)
      had_err = -1;
  }
  gt_str_array_delete(attribs);
  gt_mutex_unlock(fi->dblock);
  return had_err;
}

static int insert_feature_node(GtFeatureIndexGFFlike *fi,
                               GtFeatureNode *fn,
                               GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *toplevel;
  int had_err = 0;
  unsigned long num;

  /* TODO: locking! subgraph insertion is a transaction
      we also need to avoid nesting */

  /* collect relationships */
  gt_hashmap_reset(fi->node_to_parent_array);
  gt_feature_node_traverse_children(fn, fi->node_to_parent_array, get_parents,
                                    true, NULL);

  /* top level nodes can be pseudo nodes, check this now as
     pseudo nodes are skipped during traversal */
  if (gt_feature_node_is_pseudo(toplevel = fn)) {
    had_err = insert_single_node(fi, &num, fn, err);
  }

  /* insert nodes */
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (fn = gt_feature_node_iterator_next(fni))) {
    gt_feature_node_set_observer(fn, fi->obs);
    had_err = insert_single_node(fi, &num, fn, err);
  }

  /* nodes inserted or cached, process relationships in this subgraph */
  if (!had_err)
    gt_hashmap_foreach(fi->node_to_parent_array, set_parents, fi, err);

  /* TODO: locking! commit this subgraph to DB */

  gt_feature_node_iterator_delete(fni);
  gt_hashmap_reset(fi->node_to_parent_array);
  return had_err;
}

typedef struct {
  GtError *err;
  GtFeatureIndexGFFlike *fis;
  int had_err;
} ObserverCallbackInfo;

static void node_attribute_change_callback(GtFeatureNode *changed,
                                           GT_UNUSED bool added,
                                           GT_UNUSED const char *key,
                                           GT_UNUSED const char *value,
                                           void *data)
{
  unsigned long *id;
  ObserverCallbackInfo *info = (ObserverCallbackInfo*) data;
  gt_assert(changed);
  if ((id = node_ul_gt_hashmap_get(info->fis->cache_node2id, changed))) {
    gt_hashmap_add(info->fis->changed, changed, (void*) 1);
  }
}

static void node_attribute_delete_callback(GtFeatureNode *changed,
                                           GT_UNUSED const char *key,
                                           GT_UNUSED void *data)
{
  unsigned long *id;
  ObserverCallbackInfo *info = (ObserverCallbackInfo*) data;
  gt_assert(changed);
  if ((id = node_ul_gt_hashmap_get(info->fis->cache_node2id, changed))) {
    gt_hashmap_add(info->fis->changed, changed, (void*) 1);
  }
}

static int node_child_add_callback(GtFeatureNode *parent,
                                   GtFeatureNode *child,
                                   void *data)
{
  int had_err = 0;
  unsigned long *id, *child_id;
  ObserverCallbackInfo *info = (ObserverCallbackInfo*) data;
  gt_assert(parent && child);

  if ((id = node_ul_gt_hashmap_get(info->fis->cache_node2id, parent))) {
    if (!(child_id = node_ul_gt_hashmap_get(info->fis->cache_node2id, child))) {
      gt_hashmap_add(info->fis->added, child, (void*) 1);
    }
    gt_hashmap_add(info->fis->changed, parent, (void*) 1);
  }
  return had_err;
}

static void node_score_changed_callback(GtFeatureNode *changed,
                                        GT_UNUSED double score,
                                        void *data)
{
  unsigned long *id;
  ObserverCallbackInfo *info = (ObserverCallbackInfo*) data;
  gt_assert(changed);
  if ((id = node_ul_gt_hashmap_get(info->fis->cache_node2id, changed))) {
    gt_hashmap_add(info->fis->changed, changed, (void*) 1);
  }
}

static void node_phase_changed_callback(GtFeatureNode *changed,
                                        GT_UNUSED GtPhase phase,
                                        void *data)
{
  unsigned long *id;
  ObserverCallbackInfo *info = (ObserverCallbackInfo*) data;
  if ((id = node_ul_gt_hashmap_get(info->fis->cache_node2id, changed))) {
    gt_hashmap_add(info->fis->changed, changed, (void*) 1);
  }
}

static void node_strand_changed_callback(GtFeatureNode *changed,
                                         GT_UNUSED GtStrand strand,
                                         void *data)
{
  unsigned long *id;
  ObserverCallbackInfo *info = (ObserverCallbackInfo*) data;
  if ((id = node_ul_gt_hashmap_get(info->fis->cache_node2id, changed))) {
    gt_hashmap_add(info->fis->changed, changed, (void*) 1);
  }
}

int gt_feature_index_gfflike_add_feature_node(GtFeatureIndex *gfi,
                                              GtFeatureNode *gf,
                                              GtError *err)
{
  GtFeatureIndexGFFlike *fi;
  int had_err = 0;
  gt_assert(gfi && gf);

  fi = feature_index_gfflike_cast(gfi);
  had_err = insert_feature_node(fi,
                                (GtFeatureNode*)
                                         gt_genome_node_ref((GtGenomeNode*) gf),
                                err);
  if (!had_err)
    gt_hashmap_add(fi->ref_nodes, gf, (void*) 1);
  return had_err;
}

static int remove_node_by_id(GtFeatureIndexGFFlike *fis, unsigned long id,
                             GtError *err)
{
  unsigned long parent_count;
  int had_err = 0;

  if (had_err >= 0) {
    gt_rdb_stmt_reset(fis->stmts[GT_PSTMT_NODE_DELETE_AS_CHILD], err);
    gt_rdb_stmt_bind_int(fis->stmts[GT_PSTMT_NODE_DELETE_AS_CHILD], 0,
                         id, err);
    gt_rdb_stmt_bind_int(fis->stmts[GT_PSTMT_NODE_DELETE_AS_CHILD], 1,
                         id, err);
    had_err = gt_rdb_stmt_exec(fis->stmts[GT_PSTMT_NODE_DELETE_AS_CHILD],
                               err);
  }

  if (had_err >= 0) {
    gt_rdb_stmt_reset(fis->stmts[GT_PSTMT_GET_PARENTS_COUNT], err);
    gt_rdb_stmt_bind_int(fis->stmts[GT_PSTMT_GET_PARENTS_COUNT], 0, id, err);
    had_err = gt_rdb_stmt_exec(fis->stmts[GT_PSTMT_GET_PARENTS_COUNT], err);

    if (had_err >= 0) {
      gt_rdb_stmt_get_ulong(fis->stmts[GT_PSTMT_GET_PARENTS_COUNT],
                            0, &parent_count, err);
      if (parent_count == 0) {
        /* no more references to this node exist, we can delete it */
        /* XXX TODO: error handling */
        if (had_err >= 0) {
          gt_rdb_stmt_reset(fis->stmts[GT_PSTMT_NODE_DELETE_FEATURE],
                            err);
          gt_rdb_stmt_bind_int(fis->stmts[GT_PSTMT_NODE_DELETE_FEATURE],
                               0, id, err);
          had_err = gt_rdb_stmt_exec(fis->stmts[GT_PSTMT_NODE_DELETE_FEATURE],
                                     err);
        }
        if (had_err >= 0) {
          gt_rdb_stmt_reset(fis->stmts[GT_PSTMT_NODE_DELETE_ATTRIB_FOR_NODE],
                            err);
          gt_rdb_stmt_bind_int(fis->stmts[GT_PSTMT_NODE_DELETE_ATTRIB_FOR_NODE],
                                0, id, err);
          had_err = gt_rdb_stmt_exec(
                        fis->stmts[GT_PSTMT_NODE_DELETE_ATTRIB_FOR_NODE], err);
        }
      }
    }
  }
  if (!had_err) {
    GtFeatureNode *fn;
    if ((fn = *(GtFeatureNode**) ul_node_gt_hashmap_get(fis->cache_id2node,
                                                        id))) {
      ul_node_gt_hashmap_remove(fis->cache_id2node, id);
    }
    if (node_ul_gt_hashmap_get(fis->cache_node2id, fn)) {
      node_ul_gt_hashmap_remove(fis->cache_node2id, fn);
    }
    if (gt_hashmap_get(fis->ref_nodes, fn)) {
      gt_hashmap_remove(fis->ref_nodes, fn);
    }
  }
  return had_err;
}

int gt_feature_index_gfflike_remove_node(GtFeatureIndex *gfi,
                                         GtFeatureNode *gf,
                                         GT_UNUSED GtError *err)
{
  GtFeatureIndexGFFlike *fi;
  int had_err = 0;
  gt_assert(gfi && gf);

  fi = feature_index_gfflike_cast(gfi);
  if (gt_hashmap_get(fi->added, gf)) {
    gt_hashmap_remove(fi->added, gf);
  }
  if (gt_hashmap_get(fi->changed, gf)) {
    gt_hashmap_remove(fi->changed, gf);
  }
  gt_hashmap_add(fi->deleted, gf, (void*) 1);

  return had_err;
}

static int gt_feature_index_gfflike_save_del(void *key, GT_UNUSED void *val,
                                             void *data, GtError *err)
{
  GtFeatureNode *fn = (GtFeatureNode*) key;
  GtFeatureNodeIterator *it;
  unsigned long *id;
  ObserverCallbackInfo *oci = (ObserverCallbackInfo*) data;
  int had_err = 0;

  it = gt_feature_node_iterator_new(fn);
  id = node_ul_gt_hashmap_get(oci->fis->cache_node2id, fn);
  if (!id) {
    return 0;
  }

  while ((fn = gt_feature_node_iterator_next(it))) {
    id = node_ul_gt_hashmap_get(oci->fis->cache_node2id, fn);
    remove_node_by_id(oci->fis, *id, err);
  }
  return had_err;
}

static int gt_feature_index_gfflike_save_add(void *key, GT_UNUSED void *val,
                                             void *data, GtError *err)
{
  GtFeatureNode *fn = (GtFeatureNode*) key;
  ObserverCallbackInfo *oci = (ObserverCallbackInfo*) data;
  int had_err = 0;

  had_err = gt_feature_index_add_feature_node((GtFeatureIndex*) oci->fis, fn,
                                              err);

  return had_err;
}

typedef struct {
  GtFeatureIndexGFFlike *fis;
  unsigned long id;
  GtError *err;
} GFFlikeAttributeInfo;

static void resave_each_attribute(const char *attr_name, const char *attr_value,
                                  void *data)
{
  GFFlikeAttributeInfo *ai = (GFFlikeAttributeInfo*) data;
  gt_rdb_stmt_reset(ai->fis->stmts[GT_PSTMT_ATTRIBUTE_INSERT], ai->err);
  gt_rdb_stmt_bind_int(ai->fis->stmts[GT_PSTMT_ATTRIBUTE_INSERT], 0, ai->id,
                       ai->err);
  gt_rdb_stmt_bind_string(ai->fis->stmts[GT_PSTMT_ATTRIBUTE_INSERT], 1,
                          attr_name, ai->err);
  gt_rdb_stmt_bind_string(ai->fis->stmts[GT_PSTMT_ATTRIBUTE_INSERT], 2,
                          attr_value, ai->err);
  gt_rdb_stmt_exec(ai->fis->stmts[GT_PSTMT_ATTRIBUTE_INSERT], ai->err);
}

static int gt_feature_index_gfflike_save_chg(void *key, GT_UNUSED void *val,
                                             void *data, GtError *err)
{
  GtFeatureNode *fn = (GtFeatureNode*) key, *child;
  unsigned long *id;
  ObserverCallbackInfo *oci = (ObserverCallbackInfo*) data;
  GFFlikeAttributeInfo ai;
  GtFeatureNodeIterator *fni;
  int had_err = 0, rval;

  id = node_ul_gt_hashmap_get(oci->fis->cache_node2id, fn);
  gt_assert(id);

  /* update main feature table */
  gt_rdb_stmt_reset(oci->fis->stmts[GT_PSTMT_FEATURE_UPDATE], err);
  gt_rdb_stmt_bind_double(oci->fis->stmts[GT_PSTMT_FEATURE_UPDATE], 0,
                          gt_feature_node_score_is_defined(fn) ?
                         gt_feature_node_get_score(fn) :
                         GT_UNDEF_DOUBLE, err);
  gt_rdb_stmt_bind_int(oci->fis->stmts[GT_PSTMT_FEATURE_UPDATE], 1,
                       gt_feature_node_get_strand(fn), err);
  gt_rdb_stmt_bind_int(oci->fis->stmts[GT_PSTMT_FEATURE_UPDATE], 2,
                       gt_feature_node_get_phase(fn), err);
  gt_rdb_stmt_bind_int(oci->fis->stmts[GT_PSTMT_FEATURE_UPDATE], 3, *id,
                       err);
  rval = gt_rdb_stmt_exec(oci->fis->stmts[GT_PSTMT_FEATURE_UPDATE], err);
  if (rval < 0)
    had_err = -1;

  if (!had_err) {
    /* resave attributes: remove first... */
    gt_rdb_stmt_reset(oci->fis->stmts[GT_PSTMT_NODE_DELETE_ATTRIB_FOR_NODE],
                      err);
    gt_rdb_stmt_bind_int(oci->fis->stmts[GT_PSTMT_NODE_DELETE_ATTRIB_FOR_NODE],
                         0, *id, err);
    rval = gt_rdb_stmt_exec(
                          oci->fis->stmts[GT_PSTMT_NODE_DELETE_ATTRIB_FOR_NODE],
                          err);
    if (rval < 0)
      had_err = -1;
  }

  if (!had_err) {
    /* then insert all attributes again */
    ai.fis = oci->fis;
    ai.id = *id;
    ai.err = err;
    gt_feature_node_foreach_attribute(fn, resave_each_attribute, &ai);
  }

  if (!had_err) {
    gt_rdb_stmt_reset(oci->fis->stmts[GT_PSTMT_NODE_DELETE_AS_PARENT], err);
    gt_rdb_stmt_bind_int(oci->fis->stmts[GT_PSTMT_NODE_DELETE_AS_PARENT],
                         0, *id, err);
    rval = gt_rdb_stmt_exec(oci->fis->stmts[GT_PSTMT_NODE_DELETE_AS_PARENT],
                            err);
    fni = gt_feature_node_iterator_new_direct(fn);
    while ((child = gt_feature_node_iterator_next(fni))) {
      unsigned long *child_id;
      if ((child_id = node_ul_gt_hashmap_get(oci->fis->cache_node2id, child))) {
        gt_rdb_stmt_reset(oci->fis->stmts[GT_PSTMT_PARENT_INSERT], err);
        gt_rdb_stmt_bind_int(oci->fis->stmts[GT_PSTMT_PARENT_INSERT],
                             0, *child_id, err);
        gt_rdb_stmt_bind_int(oci->fis->stmts[GT_PSTMT_PARENT_INSERT],
                             1, *id, err);
        rval = gt_rdb_stmt_exec(oci->fis->stmts[GT_PSTMT_PARENT_INSERT],
                                err);
      }
    }
  }

  return (had_err >= 0) ? 0  : -1;
}

static int gt_feature_index_gfflike_save(GtFeatureIndex *fi,
                                         GT_UNUSED GtError *err)
{
  GtFeatureIndexGFFlike *fig;
  ObserverCallbackInfo *oci;
  int had_err = 0;
  GtRDBStmt *stmt_b, *stmt_e;
  gt_assert(fi);
  fig = feature_index_gfflike_cast(fi);
  oci = (ObserverCallbackInfo*) fig->obs->data;

  gt_mutex_lock(fig->dblock);
  stmt_b = gt_rdb_prepare(fig->db, "BEGIN TRANSACTION;", 0, err);
  stmt_e = gt_rdb_prepare(fig->db, "END TRANSACTION;", 0, err);
  gt_rdb_stmt_exec(stmt_b, err);
  if (oci && fig->deleted) {
    had_err = gt_hashmap_foreach(fig->deleted,
                                 gt_feature_index_gfflike_save_del,
                                 oci, err);
  }
  gt_hashmap_reset(fig->deleted);
  gt_rdb_stmt_exec(stmt_e, err);

  gt_rdb_stmt_reset(stmt_b, err);
  gt_rdb_stmt_exec(stmt_b, err);
  if (oci && fig->added) {
    had_err = gt_hashmap_foreach(fig->added,
                                 gt_feature_index_gfflike_save_add,
                                 oci, err);
  }
  gt_hashmap_reset(fig->added);
  gt_rdb_stmt_reset(stmt_e, err);
  gt_rdb_stmt_exec(stmt_e, err);

  gt_rdb_stmt_reset(stmt_b, err);
  gt_rdb_stmt_exec(stmt_b, err);
  if (oci && fig->changed) {
    had_err = gt_hashmap_foreach(fig->changed,
                                 gt_feature_index_gfflike_save_chg,
                                 oci, err);
  }
  gt_hashmap_reset(fig->changed);
  gt_rdb_stmt_reset(stmt_e, err);
  gt_rdb_stmt_exec(stmt_e, err);

  gt_rdb_stmt_delete(stmt_e);
  gt_rdb_stmt_delete(stmt_b);
  gt_mutex_lock(fig->dblock);

  return had_err;
}

static int get_nodes_for_stmt(GtFeatureIndexGFFlike *fi,
                              GtArray *results,
                              GtRDBStmt *stmt,
                              GtError *err)
{
  GtRDBStmt *attr_stmt, *parent_stmt;
  int had_err = 0;
  unsigned long i;
  GtArray *nodes;
  bool is_child;
  gt_assert(fi && results && stmt);
  attr_stmt = fi->stmts[GT_PSTMT_GET_ATTRIBUTE_SELECT];
  parent_stmt = fi->stmts[GT_PSTMT_GET_PARENTS_SELECT];
  nodes = gt_array_new(sizeof (unsigned long));
  HashElemInfo node_hashtype = {
    gt_ht_ptr_elem_hash,
    { NULL },
    sizeof (GtGenomeNode*),
    gt_ht_ptr_elem_cmp,
    NULL,
    NULL };
  GtHashtable *seen_as_children = gt_hashtable_new(node_hashtype);

  while (!had_err && gt_rdb_stmt_exec(stmt, err) == 0) {
    unsigned long id = GT_UNDEF_ULONG, multi_rep = GT_UNDEF_ULONG, *idp;
    int phase, is_multi = 0, strand = GT_STRAND_UNKNOWN;
    GtGenomeNode *newgn;
    GtFeatureNode *newfn;
    GtStr *seqid_str, *source_str, *type_str;
    double score;
    GtRange rng;

    seqid_str = gt_str_new();
    source_str = gt_str_new();
    type_str = gt_str_new();

    gt_rdb_stmt_get_ulong(stmt, 0, &id, err);
    gt_rdb_stmt_get_string(stmt, 1, seqid_str, err);
    gt_rdb_stmt_get_string(stmt, 2, source_str, err);
    gt_rdb_stmt_get_string(stmt, 3, type_str, err);
    gt_rdb_stmt_get_ulong(stmt, 4, &rng.start, err);
    gt_rdb_stmt_get_ulong(stmt, 5, &rng.end, err);
    gt_rdb_stmt_get_double(stmt, 6, &score, err);
    gt_rdb_stmt_get_int(stmt, 7, &strand, err);
    gt_rdb_stmt_get_int(stmt, 8, &phase, err);
    gt_rdb_stmt_get_int(stmt, 9, &is_multi, err);
    if (is_multi)
      gt_rdb_stmt_get_ulong(stmt, 10, &multi_rep, err);

    /* did we hand out or store this node already? -> reuse reference */
    if ((ul_node_gt_hashmap_get(fi->cache_id2node, id)) != NULL) {

      gt_array_add(nodes, id);

    } else {     /* otherwise build new nodes from database info */

      /* create new node */
      newgn = gt_feature_node_new(seqid_str, gt_str_get(type_str), rng.start,
                                  rng.end,strand);
      newfn = gt_feature_node_cast(newgn);
      gt_feature_node_set_phase(newfn, phase);
      gt_feature_node_set_source(newfn, source_str);
      if (score != GT_UNDEF_DOUBLE)
        gt_feature_node_set_score(newfn, score);

      /* cache node */
      ul_node_gt_hashmap_add(fi->cache_id2node, id, newfn);
      if (!(idp = node_ul_gt_hashmap_get(fi->cache_node2id, newfn)))
      {
        node_ul_gt_hashmap_add(fi->cache_node2id, newfn, id);
      }

      /* assign attributes */
      gt_rdb_stmt_reset(attr_stmt, err);
      gt_rdb_stmt_bind_ulong(attr_stmt, 0, id, err);
      while (gt_rdb_stmt_exec(attr_stmt, err) == 0) {
        GtStr *key, *value;
        key = gt_str_new();
        value = gt_str_new();
        gt_rdb_stmt_get_string(attr_stmt, 0, key, err);
        gt_rdb_stmt_get_string(attr_stmt, 1, value, err);
        gt_feature_node_set_attribute(newfn, gt_str_get(key),
                                      gt_str_get(value));
        gt_str_delete(key);
        gt_str_delete(value);
      }

      /* is this a multi-feature? */
      if (is_multi) {
        if (multi_rep == 0UL) {
          gt_feature_node_make_multi_representative(newfn);
        } else {
          GtFeatureNode *rep;
          rep = *(GtFeatureNode**) ul_node_gt_hashmap_get(fi->cache_id2node,
                                                          multi_rep);
          gt_assert(rep);
          gt_feature_node_set_multi_representative(newfn, rep);
        }
      }
      gt_array_add(nodes, id);
    }
    gt_str_delete(seqid_str);
    gt_str_delete(source_str);
    gt_str_delete(type_str);
  }

  /* rebuild DAG */
  for (i=0;i<gt_array_size(nodes);i++) {
    unsigned long id = *(unsigned long*) gt_array_get(nodes, i);
    is_child = false;
    GtFeatureNode *newfn = *(GtFeatureNode**)
                                  ul_node_gt_hashmap_get(fi->cache_id2node, id);
    gt_assert(newfn);
    /* assign parents */
    gt_rdb_stmt_reset(parent_stmt, err);
    gt_rdb_stmt_bind_ulong(parent_stmt, 0, id, err);
    while (gt_rdb_stmt_exec(parent_stmt, err) == 0) {
      unsigned long par_id;
      GtFeatureNode *parent;
      is_child = true;
      gt_rdb_stmt_get_ulong(parent_stmt, 0, &par_id, err);
      parent = *(GtFeatureNode**) ul_node_gt_hashmap_get(fi->cache_id2node,
                                                       par_id);
      gt_assert(parent);
      /* if a child has multiple parents, increase refcount */
      if (gt_hashtable_get(seen_as_children, &newfn)) {
        gt_genome_node_ref((GtGenomeNode*) newfn);
      }
      gt_feature_node_add_child(parent, newfn);
      gt_hashtable_add(seen_as_children, &newfn);
    }
    if (!is_child) {
      GtGenomeNode *newgn = (GtGenomeNode*) newfn;
      gt_array_add(results, newgn);
    }
  }
  for (i=0;i<gt_array_size(nodes);i++) {
    unsigned long id = *(unsigned long*) gt_array_get(nodes, i);
    GtFeatureNode *newfn = *(GtFeatureNode**)
                                  ul_node_gt_hashmap_get(fi->cache_id2node, id);
    gt_feature_node_set_observer(newfn, fi->obs);
  }
  gt_array_delete(nodes);
  gt_hashtable_delete(seen_as_children);
  return had_err;
}

GtArray* gt_feature_index_gfflike_get_features_for_seqid(GtFeatureIndex *gfi,
                                                         const char *seqid,
                                                         GtError *err)
{
  GtArray *a;
  GtFeatureIndexGFFlike *fi;
  GtRDBStmt *stmt;
  gt_assert(gfi && seqid);
  gt_error_check(err);
  fi = feature_index_gfflike_cast(gfi);
  stmt = fi->stmts[GT_PSTMT_GET_BY_SEQID_SELECT];
  a = gt_array_new(sizeof (GtFeatureNode*));
  gt_rdb_stmt_reset(stmt, err);
  gt_rdb_stmt_bind_string(stmt, 0, seqid, err);
  get_nodes_for_stmt(fi, a, stmt, err);
  return a;
}

int gt_feature_index_gfflike_get_features_for_range(GtFeatureIndex *gfi,
                                                    GtArray *results,
                                                    const char *seqid,
                                                    const GtRange *qry_range,
                                                    GtError *err)
{
  GtFeatureIndexGFFlike *fi;
  int retval;
  gt_error_check(err);
  GtRDBStmt *stmt;
  gt_assert(gfi && results);
  gt_error_check(err);
  fi = feature_index_gfflike_cast(gfi);
  stmt = fi->stmts[GT_PSTMT_GET_RANGE_SELECT];
  gt_mutex_lock(fi->dblock);
  gt_rdb_stmt_reset(stmt, err);
  gt_rdb_stmt_bind_string(stmt, 0, seqid, err);
  gt_rdb_stmt_bind_ulong(stmt, 1, qry_range->end, err);
  gt_rdb_stmt_bind_ulong(stmt, 2, qry_range->start, err);
  retval = get_nodes_for_stmt(fi, results, stmt, err);
  gt_mutex_unlock(fi->dblock);
  return retval;
}

int gt_feature_index_gfflike_get_all_features(GtFeatureIndex *gfi,
                                              GtArray *results,
                                              GtError *err)
{
  GtFeatureIndexGFFlike *fi;
  int retval;
  gt_error_check(err);
  GtRDBStmt *stmt;
  gt_assert(gfi && results);
  gt_error_check(err);
  fi = feature_index_gfflike_cast(gfi);
  stmt = fi->stmts[GT_PSTMT_GET_ALL];
  gt_mutex_lock(fi->dblock);
  gt_rdb_stmt_reset(stmt, err);
  retval = get_nodes_for_stmt(fi, results, stmt, err);
  gt_mutex_unlock(fi->dblock);
  return retval;
}

char* gt_feature_index_gfflike_get_first_seqid(const GtFeatureIndex *gfi,
                                               GtError *err)
{
  GtFeatureIndexGFFlike *fi;
  char *firstseqid = NULL;
  int rval, had_err = 0;
  GtStr *result;
  GtRDBStmt *stmt;
  gt_assert(gfi);
  gt_error_check(err);
  fi = feature_index_gfflike_cast((GtFeatureIndex*) gfi);
  stmt = fi->stmts[GT_PSTMT_GET_FIRST_SEQID_SELECT];
  gt_mutex_lock(fi->dblock);
  gt_rdb_stmt_reset(stmt, err);
  rval = gt_rdb_stmt_exec(stmt, err);
  result = gt_str_new();
  switch (rval) {
    case 0:
      gt_rdb_stmt_get_string(stmt, 0, result, err);
      firstseqid = gt_cstr_dup(gt_str_get(result));
      gt_assert(gt_rdb_stmt_exec(stmt, err) == 1);
      break;
    case 1:
      gt_error_set(err, "no sequence regions in index");
      had_err = -1;
      break;
    default:
      /* error is set */
      had_err = -1;
      break;
  }
  gt_str_delete(result);
  gt_mutex_unlock(fi->dblock);
  return (had_err ? NULL : firstseqid);
}

GtStrArray* gt_feature_index_gfflike_get_seqids(const GtFeatureIndex *gfi,
                                                GtError *err)
{
  GtStrArray* seqids;
  GtFeatureIndexGFFlike *fi;
  GtRDBStmt *stmt;
  int rval;
  gt_assert(gfi);
  gt_error_check(err);
  fi = feature_index_gfflike_cast((GtFeatureIndex*) gfi);
  seqids = gt_str_array_new();
  gt_mutex_lock(fi->dblock);
  stmt = fi->stmts[GT_PSTMT_GET_SEQIDS_SELECT];
  gt_rdb_stmt_reset(stmt, err);
  while ((rval = gt_rdb_stmt_exec(stmt, err)) == 0) {
    GtStr *seqid = gt_str_new();
    gt_rdb_stmt_get_string(stmt, 0, seqid, err);
    gt_str_array_add(seqids, seqid);
    gt_str_delete(seqid);
  }
  gt_mutex_unlock(fi->dblock);
  return seqids;
}

int gt_feature_index_gfflike_get_range_for_seqid(GtFeatureIndex *gfi,
                                                 GtRange *range,
                                                 const char *seqid,
                                                 GtError *err)
{
  GtFeatureIndexGFFlike *fi;
  int had_err = 0, rval;
  GtRDBStmt *stmt;
  gt_assert(gfi && range && seqid);
  gt_error_check(err);
  fi = feature_index_gfflike_cast(gfi);
  gt_mutex_lock(fi->dblock);
  stmt = fi->stmts[GT_PSTMT_GET_SEQREG_RANGE_SELECT];
  gt_rdb_stmt_reset(stmt, err);
  gt_rdb_stmt_bind_string(stmt, 0, seqid, err);
  rval = gt_rdb_stmt_exec(stmt, err);
  switch (rval) {
    case 0:
      gt_rdb_stmt_get_ulong(stmt, 0, &range->start, err);
      gt_rdb_stmt_get_ulong(stmt, 1, &range->end, err);
      gt_assert(gt_rdb_stmt_exec(stmt, err) == 1);
      break;
    case 1:
      gt_error_set(err, "sequence region '%s' does not exist", seqid);
      had_err = -1;
      break;
    default:
      /* error is set */
      had_err = -1;
      break;
  }
  gt_mutex_unlock(fi->dblock);
  return had_err;
}

int gt_feature_index_gfflike_has_seqid(const GtFeatureIndex *gfi,
                                       bool *has_seqid,
                                       const char *seqid,
                                       GtError *err)
{
  GtFeatureIndexGFFlike *fi;
  int had_err = 0;
  gt_assert(gfi);
  int rval;
  GtRDBStmt *stmt;
  gt_assert(gfi);
  gt_error_check(err);
  fi = feature_index_gfflike_cast((GtFeatureIndex*) gfi);
  gt_mutex_lock(fi->dblock);
  stmt = fi->stmts[GT_PSTMT_HAS_SEQID_SELECT];
  gt_rdb_stmt_reset(stmt, err);
  gt_rdb_stmt_bind_string(stmt, 0, seqid, err);
  rval = gt_rdb_stmt_exec(stmt, err);
  if (rval < 0) {
    /* error is set */
    had_err = -1;
  }
  if (!had_err) {
    int nof_entries;
    gt_rdb_stmt_get_int(stmt, 0, &nof_entries, err);
    *has_seqid = (nof_entries > 0);
    gt_assert(gt_rdb_stmt_exec(stmt, err) == 1);
  }
  gt_mutex_unlock(fi->dblock);
  return had_err;
}

static int
unregister_observer(void *key, GT_UNUSED void *val, GT_UNUSED void *data,
                    GT_UNUSED GtError *err)
{
  GtFeatureNode *fn = (void*) key;
  gt_assert(fn);
  gt_feature_node_unset_observer(fn);
  return 0;
}

static void unregister_observers(GtFeatureIndexGFFlike *fis, GtError *err)
{
  gt_hashmap_foreach(fis->ref_nodes, unregister_observer, NULL, err);
}

void gt_feature_index_gfflike_delete(GtFeatureIndex *gfi)
{
  GtFeatureIndexGFFlike *fi;
  unsigned long i;
  if (!gfi) return;
  fi = feature_index_gfflike_cast(gfi);
  for (i=0;i<GT_PSTMT_NOF_STATEMENTS;i++) {
    gt_rdb_stmt_delete(fi->stmts[i]);
  }
  if (fi->db)
    gt_rdb_delete(fi->db);
  gt_hashmap_delete(fi->node_to_parent_array);
  unregister_observers(fi, NULL);
  gt_feature_node_observer_delete(fi->obs);
  node_ul_gt_hashmap_delete(fi->cache_node2id);
  ul_node_gt_hashmap_delete(fi->cache_id2node);
  gt_hashmap_delete(fi->ref_nodes);
  gt_hashmap_delete(fi->seqid_cache);
  gt_hashmap_delete(fi->source_cache);
  gt_hashmap_delete(fi->string_caches);
  gt_hashmap_delete(fi->added);
  gt_hashmap_delete(fi->deleted);
  gt_hashmap_delete(fi->changed);
  gt_mutex_delete(fi->dblock);
}

const GtFeatureIndexClass* feature_index_gfflike_class(void)
{
  static const GtFeatureIndexClass *fic = NULL;
  if (!fic) {
    fic = gt_feature_index_class_new(sizeof (GtFeatureIndexGFFlike),
                                gt_feature_index_gfflike_add_region_node,
                                gt_feature_index_gfflike_add_feature_node,
                                gt_feature_index_gfflike_remove_node,
                                gt_feature_index_gfflike_get_features_for_seqid,
                                gt_feature_index_gfflike_get_features_for_range,
                                gt_feature_index_gfflike_get_first_seqid,
                                gt_feature_index_gfflike_save,
                                gt_feature_index_gfflike_get_seqids,
                                gt_feature_index_gfflike_get_range_for_seqid,
                                gt_feature_index_gfflike_has_seqid,
                                gt_feature_index_gfflike_delete);
  }
  return fic;
}

static int prepstmt_init(GtFeatureIndexGFFlike *fis, GtError *err)
{
  GtRDBStmt *r;
  gt_assert(fis);
  r = fis->stmts[GT_PSTMT_SOURCE_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT source_id FROM sources "
                        "WHERE source_name = ?",
                        1,
                        err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_SOURCE_INSERT] = gt_rdb_prepare(fis->db,
                        "INSERT INTO sources (source_name) "
                        "VALUES (?)",
                        1,
                        err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_SEQUENCEREGION_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT sequenceregion_id FROM sequenceregions "
                        "WHERE sequenceregion_name = ?",
                        1,
                        err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_SEQREG_RANGE_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT start, stop FROM sequenceregions "
                        "WHERE sequenceregion_name = ?",
                        1,
                        err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_SEQUENCEREGION_INSERT] = gt_rdb_prepare(fis->db,
                        "INSERT INTO sequenceregions "
                          "(sequenceregion_name, start, stop) "
                        "VALUES (?, ?, ?)",
                        3,
                        err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_TYPE_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT type_id FROM types "
                        "WHERE type_name = ?",
                        1,
                        err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_TYPE_INSERT] = gt_rdb_prepare(fis->db,
                        "INSERT INTO types (type_name) "
                        "VALUES (?)",
                        1,
                        err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_PARENT_INSERT] = gt_rdb_prepare(fis->db,
                        "INSERT INTO parents "
                        "VALUES (?, ?)",
                        2,
                        err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_FEATURE_INSERT] = gt_rdb_prepare(fis->db,
                        "INSERT INTO features "
                        "(seqid, source, type, start, end, score, strand, "
                        "phase, is_multi, "
                        "multi_representative, is_pseudo, is_marked) "
                        "VALUES "
                        "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                         12,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_FEATURE_UPDATE] = gt_rdb_prepare(fis->db,
                        "UPDATE features "
                        "SET score = ?, "
                        "strand = ?, "
                        "phase = ? "
                        "WHERE id = ? ",
                         4,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_ATTRIBUTE_INSERT] = gt_rdb_prepare(fis->db,
                        "INSERT INTO attributes (feature_id, keystr, value) "
                        "VALUES "
                        "(?, ?, ?)",
                         3,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_RANGE_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT f.id, s.sequenceregion_name, src.source_name, "
                        "       t.type_name, f.start, f.end, f.score, "
                        "       f.strand, f.phase, f.is_multi, "
                        "       f.multi_representative "
                        "FROM sequenceregions s, features f, "
                        "     sources src, types t "
                        "WHERE s.sequenceregion_name = ?  "
                        "AND s.sequenceregion_id = f.seqid "
                        "AND (f.start <= ? AND f.end >= ?) "
                        "AND src.source_id = f.source "
                        "AND t.type_id = f.type "
                        "ORDER BY f.id ASC",
                         3,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_ALL] = gt_rdb_prepare(fis->db,
                        "SELECT f.id, s.sequenceregion_name, src.source_name, "
                        "       t.type_name, f.start, f.end, f.score, "
                        "       f.strand, f.phase, f.is_multi, "
                        "       f.multi_representative "
                        "FROM sequenceregions s, features f, "
                        "     sources src, types t "
                        "WHERE s.sequenceregion_id = f.seqid "
                        "AND src.source_id = f.source "
                        "AND t.type_id = f.type "
                        "ORDER BY f.id ASC",
                         0,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_ATTRIBUTE_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT keystr, value FROM attributes "
                        "WHERE feature_id = ?",
                         1,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_PARENTS_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT parent FROM parents "
                        "WHERE feature_id = ?",
                         1,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_PARENTS_COUNT] = gt_rdb_prepare(fis->db,
                        "SELECT COUNT(parent) FROM parents "
                        "WHERE feature_id = ?",
                         1,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_SEQIDS_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT DISTINCT sequenceregion_name "
                        "FROM sequenceregions",
                         0,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_BY_SEQID_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT f.id, s.sequenceregion_name, src.source_name, "
                        "       t.type_name, f.start, f.end, f.score, "
                        "       f.strand, f.phase, f.is_multi, "
                        "       f.multi_representative "
                        "FROM sequenceregions s, features f, "
                        "     sources src, types t "
                        "WHERE s.sequenceregion_name = ?  "
                        "AND s.sequenceregion_id = f.seqid "
                        "AND src.source_id = f.source "
                        "AND t.type_id = f.type "
                        "ORDER BY f.id ASC",
                         1,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_HAS_SEQID_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT COUNT(sequenceregion_name) "
                        "FROM sequenceregions "
                        "WHERE sequenceregion_name = ?",
                         1,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_GET_FIRST_SEQID_SELECT] = gt_rdb_prepare(fis->db,
                        "SELECT sequenceregion_name "
                        "FROM sequenceregions "
                        "LIMIT 1",
                         0,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_DELETE_FEATURE] = gt_rdb_prepare(fis->db,
                        "DELETE FROM features "
                        "WHERE id = ? ",
                         1,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_DELETE_AS_CHILD] = gt_rdb_prepare(fis->db,
                        "DELETE FROM parents "
                        "WHERE feature_id = ? "
                        "OR parent = ?",
                         2,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_DELETE_AS_PARENT] = gt_rdb_prepare(fis->db,
                        "DELETE FROM parents "
                        "WHERE parent = ?",
                         1,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_UPDATE_SCORE] = gt_rdb_prepare(fis->db,
                        "UPDATE features "
                        "SET score = ? "
                        "WHERE id = ? ",
                         2,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_UPDATE_STRAND] = gt_rdb_prepare(fis->db,
                        "UPDATE features "
                        "SET strand = ? "
                        "WHERE id = ? ",
                         2,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_UPDATE_PHASE] = gt_rdb_prepare(fis->db,
                        "UPDATE features "
                        "SET phase = ? "
                        "WHERE id = ? ",
                         2,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_UPDATE_ATTRIB] = gt_rdb_prepare(fis->db,
                        "UPDATE attributes "
                        "SET value = ? "
                        "WHERE feature_id = ? "
                        "AND keystr = ?",
                         3,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_DELETE_ATTRIB] = gt_rdb_prepare(fis->db,
                        "DELETE FROM attributes "
                        "WHERE feature_id = ? "
                        "AND keystr = ?",
                         2,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_DELETE_ATTRIB_FOR_NODE] = gt_rdb_prepare(fis->db,
                        "DELETE FROM attributes "
                        "WHERE feature_id = ? ",
                         1,
                         err);
  if (!r) return -1;
  r = fis->stmts[GT_PSTMT_NODE_ADD_CHILD] = gt_rdb_prepare(fis->db,
                        "INSERT INTO parents (feature_id, parent) "
                        "VALUES (?, ?)",
                         2,
                         err);
  if (!r) return -1;
  return 0;
}

static void delete_ref_node(GtGenomeNode *node)
{
  if (!node) return;
  gt_genome_node_delete(node);
}

GtFeatureIndex* anno_db_gfflike_build(GtAnnoDBSchema *schema, GtRDB *db,
                                      GtError *err)
{
  int had_err = 0;
  GtFeatureIndex *fi = NULL;
  GtFeatureIndexGFFlike *fis;
  GtAnnoDBGFFlike *adg;
  ObserverCallbackInfo *oci;
  gt_assert(schema && db);
  gt_error_check(err);

  adg = anno_db_gfflike_cast(schema);
  had_err = gt_rdb_accept(db, adg->visitor, err);

  if (!had_err) {
    fi = gt_feature_index_create(feature_index_gfflike_class());
    fis = feature_index_gfflike_cast(fi);
    fis->obs = gt_feature_node_observer_new();
    fis->obs->data = oci = gt_calloc(1, sizeof (ObserverCallbackInfo));
    fis->node_to_parent_array = gt_hashmap_new(GT_HASH_DIRECT, NULL,
                                               (GtFree) gt_array_delete);
    fis->cache_node2id = node_ul_gt_hashmap_new();
    fis->cache_id2node = ul_node_gt_hashmap_new();
    fis->string_caches = gt_hashmap_new(GT_HASH_DIRECT,
                                        NULL,
                                        (GtFree) gt_hashmap_delete);
    fis->seqid_cache = gt_hashmap_new(GT_HASH_STRING, NULL,
                                      (GtFree) gt_str_delete);
    fis->ref_nodes   = gt_hashmap_new(GT_HASH_DIRECT, (GtFree) delete_ref_node,
                                      NULL);
    fis->source_cache = gt_hashmap_new(GT_HASH_STRING, NULL,
                                       (GtFree) gt_str_delete);
    fis->dblock = gt_mutex_new();

    /* set up callbacks */
    oci->fis = fis;
    fis->changed = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    fis->deleted = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    fis->added = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    oci->had_err = 0;
    oci->err = err;
    fis->obs->deleted = NULL;
    fis->obs->score_changed = node_score_changed_callback;
    fis->obs->strand_changed = node_strand_changed_callback;
    fis->obs->phase_changed = node_phase_changed_callback;
    fis->obs->attribute_changed = node_attribute_change_callback;
    fis->obs->attribute_deleted = node_attribute_delete_callback;
    fis->obs->child_added = node_child_add_callback;
    fis->db = gt_rdb_ref(db);

    if (prepstmt_init(fis, err)) {
      gt_feature_index_delete(fi);
      fi = NULL;
    }
  }
  return fi;
}

const GtAnnoDBSchemaClass* gt_anno_db_gfflike_class()
{
  static const GtAnnoDBSchemaClass *adbsc = NULL;
  gt_class_alloc_lock_enter();
  if (!adbsc) {
    adbsc = gt_anno_db_schema_class_new(sizeof (GtAnnoDBGFFlike),
                                        anno_db_gfflike_free,
                                        anno_db_gfflike_build);
  }
  gt_class_alloc_lock_leave();
  return adbsc;
}

static const GtRDBVisitorClass* gfflike_setup_visitor_class()
{
  static const GtRDBVisitorClass *svc = NULL;
  gt_class_alloc_lock_enter();
  if (!svc) {
    svc = gt_rdb_visitor_class_new(sizeof (GtAnnoDBGFFlike),
                                   NULL,
                                   anno_db_gfflike_init_sqlite,
                                   anno_db_gfflike_init_mysql);
  }
  gt_class_alloc_lock_leave();
  return svc;
}

static GtRDBVisitor* gfflike_setup_visitor_new(GtAnnoDBGFFlike *adb)
{
  GtRDBVisitor *v = gt_rdb_visitor_create(gfflike_setup_visitor_class());
  GFFlikeSetupVisitor *sv = gfflike_setup_visitor_cast(v);
  gt_assert(adb);
  sv->annodb = adb;
  return v;
}

GtAnnoDBSchema* gt_anno_db_gfflike_new(void)
{
  GtAnnoDBSchema *s = gt_anno_db_schema_create(gt_anno_db_gfflike_class());
  GtAnnoDBGFFlike *adg = anno_db_gfflike_cast(s);
  adg->visitor = gfflike_setup_visitor_new(adg);
  return s;
}

int gt_anno_db_gfflike_unit_test(GtError *err)
{
  int had_err = 0, status = 0;
  GtFeatureIndex *fi = NULL;
  GtError *testerr = NULL;
  GtAnnoDBSchema *adb = NULL;
  FILE *tmpfp;
#ifdef HAVE_SQLITE
  GtRDB *rdb;
#endif
  GtStr* tmpfilename;
  gt_error_check(err);

  testerr = gt_error_new();

  tmpfilename = gt_str_new();
  tmpfp = gt_xtmpfp(tmpfilename);
  gt_fa_xfclose(tmpfp);

#ifdef HAVE_SQLITE
  rdb = gt_rdb_sqlite_new(gt_str_get(tmpfilename), testerr);
  gt_ensure(had_err, rdb != NULL);
  if (!had_err) {
    adb = gt_anno_db_gfflike_new();
    gt_ensure(had_err, adb != NULL);
  }

  if (!had_err) {
    fi = gt_anno_db_schema_get_feature_index(adb, rdb, testerr);
    gt_ensure(had_err, fi != NULL);
  }
#endif

  if (!had_err) {
    /* run generic feature index tests */
    status = gt_feature_index_unit_test(fi, testerr);
    gt_ensure(had_err, status == 0);
  }

  gt_xremove(gt_str_get(tmpfilename));
  gt_str_delete(tmpfilename);
  gt_feature_index_delete(fi);
  gt_anno_db_schema_delete(adb);
#ifdef HAVE_SQLITE
  gt_rdb_delete((GtRDB*) rdb);
#endif
  gt_error_delete(testerr);
  return had_err;
}
