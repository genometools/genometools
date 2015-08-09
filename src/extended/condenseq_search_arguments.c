/*
  Copyright (c) 2015 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include "core/logger_api.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/str_api.h"
#include "extended/condenseq.h"

#include "extended/condenseq_search_arguments.h"

struct GtCondenseqSearchArguments {
  GtStr *dbpath;
  bool   verbose;
};

GtCondenseqSearchArguments *gt_condenseq_search_arguments_new(void)
{
  GtCondenseqSearchArguments *csi = gt_malloc(sizeof(*csi));
  csi->dbpath = gt_str_new();
  return csi;
}

void gt_condenseq_search_arguments_delete(
                         GtCondenseqSearchArguments *condenseq_search_arguments)
{
  if (condenseq_search_arguments != NULL) {
    gt_str_delete(condenseq_search_arguments->dbpath);
    gt_free(condenseq_search_arguments);
  }
}

GtCondenseq *gt_condenseq_search_arguments_read_condenseq(
                   const GtCondenseqSearchArguments *condenseq_search_arguments,
                   GtLogger *logger,
                   GtError *err)
{
  return gt_condenseq_new_from_file(
                                 gt_str_get(condenseq_search_arguments->dbpath),
                                 logger, err);
}

bool gt_condenseq_search_arguments_verbose(
                   const GtCondenseqSearchArguments *condenseq_search_arguments)
{
  return condenseq_search_arguments->verbose;
}

void gt_condenseq_search_register_options(
                         GtCondenseqSearchArguments *condenseq_search_arguments,
                         GtOptionParser *option_parser)
{
  GtOption *option;
   /* -db */
  option = gt_option_new_filename("db", "path of (compressed) fasta database",
                                  condenseq_search_arguments->dbpath);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(option_parser, option);

  /* -verbose */
  option = gt_option_new_bool("verbose", "verbose output",
                              &condenseq_search_arguments->verbose, false);
  gt_option_parser_add_option(option_parser, option);
}

GtStr *gt_condenseq_search_arguments_db_filename(
                         GtCondenseqSearchArguments *condenseq_search_arguments,
                         const char *suffix)
{
  GtStr *file_path = gt_str_clone(condenseq_search_arguments->dbpath);
  gt_str_append_cstr(file_path, suffix);
  return file_path;
}
