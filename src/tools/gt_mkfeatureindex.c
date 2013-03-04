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
#include "core/str_array_api.h"
#include "core/unused_api.h"
#include "core/xposix.h"
#include "extended/anno_db_gfflike_api.h"
#include "extended/bed_in_stream.h"
#include "extended/feature_index_api.h"
#include "extended/feature_stream_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gtf_in_stream.h"
#include "extended/rdb_api.h"
#ifdef HAVE_MYSQL
#include "core/password_entry.h"
#include "extended/rdb_mysql_api.h"
#endif
#include "extended/rdb_sqlite_api.h"
#include "tools/gt_mkfeatureindex.h"

#define GT_SQLITE_BACKEND_STRING "sqlite"
#define GT_MYSQL_BACKEND_STRING  "mysql"

typedef struct {
  GtStr *backend,
        *filename,
        *host,
        *user,
        *pass,
        *database,
        *input;
  int port;
  bool verbose,
       force;
} GtMkfeatureindexArguments;

static void* gt_mkfeatureindex_arguments_new(void)
{
  GtMkfeatureindexArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->backend = gt_str_new();
  arguments->filename = gt_str_new();
  arguments->host = gt_str_new();
  arguments->user = gt_str_new();
  arguments->pass = gt_str_new();
  arguments->database = gt_str_new();
  arguments->input = gt_str_new();
  return arguments;
}

static void gt_mkfeatureindex_arguments_delete(void *tool_arguments)
{
  GtMkfeatureindexArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->backend);
  gt_str_delete(arguments->filename);
  gt_str_delete(arguments->host);
  gt_str_delete(arguments->user);
  gt_str_delete(arguments->pass);
  gt_str_delete(arguments->database);
  gt_str_delete(arguments->input);
  gt_free(arguments);
}

static GtOptionParser* gt_mkfeatureindex_option_parser_new(void *tool_arguments)
{
  GtMkfeatureindexArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *backend_option, *filenameoption;
  static const char *backends[] = {
    GT_SQLITE_BACKEND_STRING,
#ifdef HAVE_MYSQL
    GT_MYSQL_BACKEND_STRING,
#endif
    NULL
  };
  static const char *inputs[] = {
    "gff",
    "bed",
    "gtf",
    NULL
  };
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] <input filename(s)>",
                            "Creates a new FeatureIndex from annotation data.");

  /* -force */
  option = gt_option_new_bool("force", "force writing to output file",
                              &arguments->force, false);
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

  /* -input */
  option = gt_option_new_choice("input", "input data format\n"
                                       "choose from gff|bed|gtf",
                             arguments->input, inputs[0], inputs);
  gt_option_parser_add_option(op, option);

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

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_mkfeatureindex_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GT_UNUSED GtMkfeatureindexArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

static int gt_mkfeatureindex_runner(int argc,
                                    const char **argv,
                                    int parsed_args,
                                    void *tool_arguments,
                                    GtError *err)
{
  GtMkfeatureindexArguments *arguments = tool_arguments;
  GtNodeStream *in_stream = NULL,
              *feature_stream = NULL;
  GtRDB *rdb = NULL;
  GtAnnoDBSchema *adb = NULL;
  GtFeatureIndex *fis = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

#ifdef HAVE_SQLITE
  if (strcmp(gt_str_get(arguments->backend),
             GT_SQLITE_BACKEND_STRING) == 0) {
    if (gt_file_exists(gt_str_get(arguments->filename))) {
      if (arguments->force) {
        gt_xunlink(gt_str_get(arguments->filename));
      } else {
        gt_error_set(err, "file \"%s\" exists already. use option -force to "
                     "overwrite", gt_str_get(arguments->filename));
        had_err = -1;
      }
    }
    if (!had_err) {
      rdb = gt_rdb_sqlite_new(gt_str_get(arguments->filename), err);
      if (!rdb)
        had_err = -1;
    }
  }
#endif
#ifdef HAVE_MYSQL
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
#endif

  adb = gt_anno_db_gfflike_new();
  if (!had_err && !adb)
    had_err = -1;

  if (!had_err) {
    fis = gt_anno_db_schema_get_feature_index(adb, rdb, err);
    if (!fis)
      had_err = -1;
  }

  if (!had_err) {
    if (strcmp(gt_str_get(arguments->input), "gff") == 0)
    {
      in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                 argv + parsed_args);
      if (arguments->verbose)
        gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) in_stream);
    } else if (strcmp(gt_str_get(arguments->input), "bed") == 0)
    {
      if (argc - parsed_args == 0)
        in_stream = gt_bed_in_stream_new(NULL);
      else
        in_stream = gt_bed_in_stream_new(argv[parsed_args]);
    } else if (strcmp(gt_str_get(arguments->input), "gtf") == 0)
    {
      if (argc - parsed_args == 0)
        in_stream = gt_gtf_in_stream_new(NULL);
      else
        in_stream = gt_gtf_in_stream_new(argv[parsed_args]);
    }
    gt_assert(in_stream);

    feature_stream = gt_feature_stream_new(in_stream, fis);
    had_err = gt_node_stream_pull(feature_stream, err);
  }
  gt_node_stream_delete(feature_stream);
  gt_node_stream_delete(in_stream);
  gt_feature_index_delete(fis);
  gt_anno_db_schema_delete(adb);
  gt_rdb_delete(rdb);
  return had_err;
}

GtTool* gt_mkfeatureindex(void)
{
  return gt_tool_new(gt_mkfeatureindex_arguments_new,
                     gt_mkfeatureindex_arguments_delete,
                     gt_mkfeatureindex_option_parser_new,
                     gt_mkfeatureindex_arguments_check,
                     gt_mkfeatureindex_runner);
}
