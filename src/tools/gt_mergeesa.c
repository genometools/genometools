/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "core/error.h"
#include "core/logger.h"
#include "core/option_api.h"
#include "core/versionfunc.h"
#include "match/test-mergeesa.pr"
#include "tools/gt_mergeesa.h"

static GtOPrval parse_options(GtStr *indexname,GtStrArray *indexnametab,
                              int *parsed_args, int argc,
                              const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  GtOption *option;

  gt_error_check(err);
  op = gt_option_parser_new("storeindex <mkvindex1> <mkvindex2> ...",
                         "Merge indexes into one index.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");
  option = gt_option_new_filename_array("ii",
                                    "specify input index files (mandatory)",
                                    indexnametab);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("indexname",
                             "specify index to be created",
                             indexname, NULL);

  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_mergeesa(int argc, const char **argv, GtError *err)
{
  GtStr *storeindex;
  GtStrArray *indexnametab;
  bool haserr = false;
  int parsed_args;

  gt_error_check(err);

  storeindex = gt_str_new();
  indexnametab = gt_str_array_new();
  switch (parse_options(storeindex, indexnametab, &parsed_args, argc, argv,
                        err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
         haserr = true; break;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }
  if (!haserr)
  {
    unsigned long i;
    GtLogger *logger;

    printf("# storeindex=%s\n",gt_str_get(storeindex));
    for (i=0; i<gt_str_array_size(indexnametab); i++)
    {
      printf("# input=%s\n",gt_str_array_get(indexnametab,i));
    }
    logger = gt_logger_new(false, GT_LOGGER_DEFLT_PREFIX, stdout);
    if (gt_performtheindexmerging(storeindex,
                              indexnametab,
                              logger,
                              err) != 0)
    {
      haserr = true;
    }
    gt_logger_delete(logger);
  }
  gt_str_delete(storeindex);
  gt_str_array_delete(indexnametab);
  return haserr ? -1 : 0;
}
