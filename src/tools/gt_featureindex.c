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

#include <string.h>
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/password_entry.h"
#include "core/str_api.h"
#include "core/unused_api.h"
#include "extended/anno_db_gfflike_api.h"
#include "extended/anno_db_schema_api.h"
#include "extended/feature_index_api.h"
#include "extended/feature_node.h"
#include "extended/feature_stream_api.h"
#include "extended/gff3_visitor.h"
#include "extended/rdb_api.h"
#ifdef HAVE_MYSQL
#include "extended/rdb_mysql_api.h"
#endif
#ifdef HAVE_SQLITE
#include "extended/rdb_sqlite_api.h"
#endif
#include "tools/gt_featureindex.h"

#define GT_SQLITE_BACKEND_STRING "sqlite"
#define GT_MYSQL_BACKEND_STRING  "mysql"

typedef struct {
  GtRange qry_rng;
  GtStr *seqid;
  GtStr *backend,
        *filename,
        *host,
        *user,
        *pass,
        *database;
  int port;
  bool verbose,
       retain,
       child_callback_check,
       attributes_callback_check;
  GtOption *rngopt;
} GtFeatureindexArguments;

static void* gt_featureindex_arguments_new(void)
{
  GtFeatureindexArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->seqid = gt_str_new();
  arguments->backend = gt_str_new();
  arguments->filename = gt_str_new();
  arguments->host = gt_str_new();
  arguments->user = gt_str_new();
  arguments->database = gt_str_new();
  return arguments;
}

static void gt_featureindex_arguments_delete(void *tool_arguments)
{
  GtFeatureindexArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->seqid);
  gt_str_delete(arguments->backend);
  gt_str_delete(arguments->filename);
  gt_str_delete(arguments->host);
  gt_str_delete(arguments->user);
  gt_str_delete(arguments->pass);
  gt_str_delete(arguments->database);
  gt_option_delete(arguments->rngopt);
  gt_free(arguments);
}

static GtOptionParser* gt_featureindex_option_parser_new(void *tool_arguments)
{
  GtFeatureindexArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *backend_option, *filenameoption;
  static const char *backends[] = {
#ifdef HAVE_SQLITE
    GT_SQLITE_BACKEND_STRING,
#endif
#ifdef HAVE_MYSQL
    GT_MYSQL_BACKEND_STRING,
#endif
    NULL
  };
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] indexfilename",
                            "Retrieve annotations from a persistent "
                            "feature index as GFF3 output.");

  arguments->rngopt = gt_option_new_range("range",
                                          "range constraint for index query",
                                          &arguments->qry_rng, NULL);
  gt_option_parser_add_option(op, gt_option_ref(arguments->rngopt));

  option = gt_option_new_string("seqid", "sequence region",
                                arguments->seqid, NULL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("retain", "retain attributes",
                              &arguments->retain, true);
  gt_option_parser_add_option(op, option);

  /* -backend */
  backend_option = gt_option_new_choice("backend", "database backend to use\n"
                                        "choose from ["
#ifdef HAVE_SQLITE
                                         GT_SQLITE_BACKEND_STRING
#endif
#ifdef HAVE_MYSQL
                                        "|" GT_MYSQL_BACKEND_STRING
#endif
                                        "]",
                                        arguments->backend, backends[0],
                                        backends);
  gt_option_parser_add_option(op, backend_option);

  /* -filename */
  filenameoption = gt_option_new_string("filename",
                                        "filename for feature database "
                                        "(sqlite backend only)",
                                        arguments->filename, NULL);
  gt_option_parser_add_option(op, filenameoption);

#ifdef HAVE_MYSQL
  /* -host */
  option = gt_option_new_string("host", "hostname for database connection",
                                arguments->host, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory_either(option, filenameoption);

  /* -user */
  option = gt_option_new_string("user", "username for database connection",
                                arguments->user, NULL);
  gt_option_parser_add_option(op, option);

  /* -database */
  option = gt_option_new_string("database",
                                "database name for database connection",
                                arguments->database, NULL);
  gt_option_parser_add_option(op, option);

  /* -port */
  option = gt_option_new_int_max("port", "port for database connection",
                                 &arguments->port, 3333, 65534);
  gt_option_parser_add_option(op, option);
#else
  gt_option_is_mandatory(filenameoption);
#endif

  option = gt_option_new_bool("child_check", "test callbacks for "
                                             "child node addition",
                              &arguments->child_callback_check, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  option = gt_option_new_bool("attributes_check", "test callbacks for "
                                                  "attribute modification",
                              &arguments->attributes_callback_check, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_featureindex_arguments_check(GT_UNUSED int rest_argc,
                                           GT_UNUSED void *tool_arguments,
                                           GT_UNUSED GtError *err)
{
  GT_UNUSED GtFeatureindexArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

static int gt_featureindex_runner(GT_UNUSED int argc,
                                  GT_UNUSED const char **argv,
                                  GT_UNUSED int parsed_args,
                                  void *tool_arguments,
                                  GtError *err)
{
  GtFeatureindexArguments *arguments = tool_arguments;
  GtFeatureIndex *fi = NULL;
  GtArray *results = NULL;
  GtRange rng;
  GtRDB *rdb = NULL;
  GtAnnoDBSchema *adbs = NULL;
  GtNodeVisitor *gff3visitor = NULL;
  GtGenomeNode *regn = NULL;
  unsigned long i = 0;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

#ifdef HAVE_SQLITE
  if (!had_err) {
    if (strcmp(gt_str_get(arguments->backend),
               GT_SQLITE_BACKEND_STRING) == 0) {
      if (!gt_file_exists(gt_str_get(arguments->filename))) {
        gt_error_set(err, "file '%s' does not exist",
                     gt_str_get(arguments->filename));
        had_err = -1;
      }
      if (!had_err) {
        rdb = gt_rdb_sqlite_new(gt_str_get(arguments->filename), err);
        if (!rdb)
          had_err = -1;
      }
    }
  }
#endif
#ifdef HAVE_MYSQL
  if (!had_err) {
    if (strcmp(gt_str_get(arguments->backend),
                      GT_MYSQL_BACKEND_STRING) == 0) {
      arguments->pass = gt_get_password("password: ", err);
      rdb = gt_rdb_mysql_new(gt_str_get(arguments->host),
                             arguments->port,
                             gt_str_get(arguments->database),
                             gt_str_get(arguments->user),
                             gt_str_get(arguments->pass),
                             err);
      if (!rdb)
        had_err = -1;
    }
  }
#endif
  if (!had_err)
    adbs = gt_anno_db_gfflike_new();

  if (!had_err && !adbs)
    had_err = -1;

  if (!had_err) {
    fi = gt_anno_db_schema_get_feature_index(adbs, rdb, err);
    had_err = fi ? 0 : -1;
  }

  if (!had_err && gt_str_length(arguments->seqid) == 0) {
    char *firstseqid = gt_feature_index_get_first_seqid(fi, err);
    if (firstseqid == NULL)
      had_err = -1;
    else {
      gt_str_append_cstr(arguments->seqid, firstseqid);
      gt_free(firstseqid);
    }
  }

  if (!had_err && !gt_option_is_set(arguments->rngopt)) {
    had_err = gt_feature_index_get_range_for_seqid(fi, &arguments->qry_rng,
                                                   gt_str_get(arguments->seqid),
                                                   err);
  }

  if (!had_err) {
    results = gt_array_new(sizeof (GtFeatureNode*));
    had_err = gt_feature_index_get_features_for_range(fi, results,
                                                   gt_str_get(arguments->seqid),
                                                   &arguments->qry_rng, err);
  }
  if (!had_err) {
    gff3visitor = gt_gff3_visitor_new(NULL);
    if (arguments->retain)
      gt_gff3_visitor_retain_id_attributes((GtGFF3Visitor*) gff3visitor);

    had_err = gt_feature_index_get_range_for_seqid(fi, &rng,
                                                   gt_str_get(arguments->seqid),
                                                   err);
  }
  if (!had_err) {
    regn = gt_region_node_new(arguments->seqid, rng.start, rng.end);
    gt_genome_node_accept(regn, gff3visitor, err);
    gt_genome_node_delete(regn);

    for (i=0; i<gt_array_size(results); i++) {
      GtGenomeNode *gn = *(GtGenomeNode**) gt_array_get(results, i);
      GtRange nrng = gt_genome_node_get_range(gn);

      if (gt_feature_node_try_cast(gn)) {
        GtFeatureNode *fgn = (GtFeatureNode*) gn;
        if (arguments->child_callback_check) {
          GtFeatureNode *fn;
          fn = (GtFeatureNode*) gt_feature_node_new(arguments->seqid, "test",
                                                    nrng.start, nrng.end,
                                                    GT_STRAND_FORWARD);
          gt_feature_node_add_child((GtFeatureNode*) gn, fn);
        }
        if (arguments->attributes_callback_check) {
          if (gt_feature_node_get_attribute(fgn, "foo")) {
            gt_feature_node_set_attribute(fgn, "foo", "boooo");
            gt_feature_node_remove_attribute(fgn, "foo");
          }
          gt_feature_node_add_attribute(fgn, "foo", "bar");
        }
      }
      gt_genome_node_accept(gn, gff3visitor, err);
      gt_genome_node_delete(gn);
    }
  }

  gt_array_delete(results);
  gt_node_visitor_delete(gff3visitor);
  gt_rdb_delete(rdb);
  gt_anno_db_schema_delete(adbs);

  gt_feature_index_delete(fi);
  return had_err;
}

GtTool* gt_featureindex(void)
{
  return gt_tool_new(gt_featureindex_arguments_new,
                  gt_featureindex_arguments_delete,
                  gt_featureindex_option_parser_new,
                  gt_featureindex_arguments_check,
                  gt_featureindex_runner);
}
