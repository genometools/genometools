/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#include <stdio.h>
#include <string.h>
#include "core/error.h"
#include "core/logger.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/versionfunc.h"
#include "match/eis-bwtseq.h"
#include "core/encseq.h"
#include "match/sarr-def.h"
#include "match/eis-bwtseq-param.h"
#include "tools/gt_packedindex_trsuftab.h"

#define DEFAULT_PROGRESS_INTERVAL  100000UL

struct trSufTabOptions
{
  struct bwtOptions idx;
  bool verboseOutput;
};

static GtOPrval
parseTrSufTabOptions(int *parsed_args, int argc, const char **argv,
                     struct trSufTabOptions *params, const GtStr *projectName,
                     GtError *err);

extern int
gt_packedindex_trsuftab(int argc, const char *argv[], GtError *err)
{
  struct trSufTabOptions params;
  BWTSeq *bwtSeq = NULL;
  GtStr *inputProject = NULL;
  int parsedArgs;
  bool had_err = false;
  GtLogger *logger = NULL;
  inputProject = gt_str_new();

  do {
    gt_error_check(err);
    {
      bool exitNow = false;
      switch (parseTrSufTabOptions(&parsedArgs, argc, argv, &params,
                                   inputProject, err))
      {
      case GT_OPTION_PARSER_OK:
        break;
      case GT_OPTION_PARSER_ERROR:
        had_err = true;
        exitNow = true;
        break;
      case GT_OPTION_PARSER_REQUESTS_EXIT:
        exitNow = true;
        break;
      }
      if (exitNow)
        break;
    }
    gt_str_set(inputProject, argv[parsedArgs]);
    logger = gt_logger_new(params.verboseOutput,
                           GT_LOGGER_DEFLT_PREFIX, stdout);
    bwtSeq = gt_trSuftab2BWTSeq(&params.idx.final, logger, err);
    had_err = bwtSeq == NULL;
    if (had_err)
      break;
  } while (0);
  if (bwtSeq) gt_deleteBWTSeq(bwtSeq);
  if (logger) gt_logger_delete(logger);
  if (inputProject) gt_str_delete(inputProject);
  return had_err?-1:0;
}

static GtOPrval
parseTrSufTabOptions(int *parsed_args, int argc, const char **argv,
                     struct trSufTabOptions *params, const GtStr *projectName,
                     GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  GtOption *option;

  gt_error_check(err);
  op = gt_option_parser_new("indexname",
                         "Build BWT packedindex for project <indexname>.");

  gt_registerPackedIndexOptions(op, &params->idx, BWTDEFOPT_MULTI_QUERY,
                             projectName);

  option = gt_option_new_bool("v",
                           "print verbose progress information",
                           &params->verboseOutput,
                           false);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 1, 1);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  /* compute parameters currently not set from command-line or
   * determined indirectly */
  gt_computePackedIndexDefaults(&params->idx, BWTBaseFeatures);

  gt_option_parser_delete(op);

  return oprval;
}
