/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "core/error_api.h"
#include "core/str_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/option_api.h"
#include "core/tool_api.h"
#include "core/versionfunc.h"
#include "match/esa-seqread.h"
#include "match/esa-shulen.h"
#include "match/esa-map.h"
#include "tools/gt_shulen.h"

typedef struct
{
  bool scanfile, beverbose;
  GtStr *indexname;
  GtStrArray *queryfilenames;
} Shulengthoptions;

static int callmultishulengthdist(const char *indexname,
                                  bool scanfile,
                                  GtLogger *logger,
                                  GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(indexname,
                                                   SARR_LCPTAB |
                                                   SARR_SUFTAB |
                                                   SARR_ESQTAB,
                                                   scanfile
                                                    ? SEQ_scan : SEQ_mappedboth,
                                                   logger,
                                                   err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_multiesa2shulengthdist_print(ssar,
                                  gt_encseqSequentialsuffixarrayreader(ssar),
                                  err) != 0)
    {
      haserr = true;
    }
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}

static int callpairswisesshulendistdist(const char *indexname,
                                        const GtStrArray *queryfilenames,
                                        GtLogger *logger,
                                        GtError *err)
{
  bool haserr = false;
  Suffixarray suffixarray;

  if (gt_mapsuffixarray(&suffixarray,
                        SARR_SUFTAB |
                        SARR_ESQTAB,
                        indexname,
                        logger,
                        err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    unsigned long totalgmatchlength = 0;

    if (gt_esa2shulengthqueryfiles(&totalgmatchlength,
                                   &suffixarray,
                                   queryfilenames,
                                   err) != 0)
    {
      haserr = true;
    } else
    {
      printf("%lu\n",totalgmatchlength);
    }
  }
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static void *gt_shulengthdist_arguments_new(void)
{
  Shulengthoptions *arguments;

  arguments = gt_malloc(sizeof (*arguments));
  arguments->indexname = gt_str_new();
  arguments->queryfilenames = gt_str_array_new();
  return arguments;
}

static void gt_shulengthdist_arguments_delete(void *tool_arguments)
{
  Shulengthoptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->indexname);
  gt_str_array_delete(arguments->queryfilenames);
  gt_free(arguments);
}

static GtOptionParser *gt_shulengthdist_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option, *queryoption, *scanoption;
  Shulengthoptions *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -ii indexname",
                            "Compute distribution of pairwise "
                            "shustring lengths.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  option = gt_option_new_string("ii",
                                "Specify input index",
                                arguments->indexname, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  scanoption = gt_option_new_bool("scan","scan index rather than mapping "
                                  "it to main memory",
                                  &arguments->scanfile,
                                  false);
  gt_option_parser_add_option(op, scanoption);

  queryoption = gt_option_new_filename_array("q",
                                             "Specify query files",
                                             arguments->queryfilenames);
  gt_option_is_development_option(queryoption);
  gt_option_parser_add_option(op, queryoption);

  option = gt_option_new_bool("v",
                              "be verbose ",
                              &arguments->beverbose,
                              false);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(queryoption,scanoption);
  return op;
}

static int gt_shulengthdist_arguments_check(GT_UNUSED int rest_argc,
                                            GT_UNUSED void *tool_arguments,
                                            GT_UNUSED GtError *err)
{
  return 0;
}

static int gt_shulengthdist_runner(GT_UNUSED int argc,
                                   GT_UNUSED const char **argv,
                                   GT_UNUSED int parsed_args,
                                   void *tool_arguments, GtError *err)
{
  bool haserr = false;
  Shulengthoptions *arguments = tool_arguments;
  GtLogger *logger = NULL;

  gt_error_check(err);
  logger = gt_logger_new(arguments->beverbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (parsed_args < argc)
  {
    gt_error_set(err,"superfluous arguments: \"%s\"",argv[argc-1]);
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_str_array_size(arguments->queryfilenames) == 0)
    {
      if (callmultishulengthdist(gt_str_get(arguments->indexname),
                                 arguments->scanfile,
                                 logger,
                                 err) != 0)
      {
        haserr = true;
      }
    } else
    {
      if (callpairswisesshulendistdist(gt_str_get(arguments->indexname),
                                       arguments->queryfilenames,
                                       logger,
                                       err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_logger_delete(logger);
  return haserr ? -1 : 0;
}

GtTool* gt_shulengthdist(void)
{
  return gt_tool_new(gt_shulengthdist_arguments_new,
                     gt_shulengthdist_arguments_delete,
                     gt_shulengthdist_option_parser_new,
                     gt_shulengthdist_arguments_check,
                     gt_shulengthdist_runner);
}
