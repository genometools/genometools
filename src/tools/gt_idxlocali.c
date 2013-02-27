/*
  Copyright (c) 2009 Z. Tang <tangzhihao0117@hotmail.com>
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include "core/option_api.h"
#include "core/ma.h"
#include "core/str_array.h"
#include "core/unused_api.h"
#include "core/tool_api.h"
#include "match/idxlocali.h"
#include "tools/gt_idxlocali.h"

static void *gt_idxlocali_arguments_new(void)
{
  return gt_malloc(sizeof (IdxlocaliOptions));
}

static void gt_idxlocali_arguments_delete(void *tool_arguments)
{
  IdxlocaliOptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->indexname);
  gt_str_array_delete(arguments->queryfiles);
  gt_option_delete(arguments->refoptionesaindex);
  gt_option_delete(arguments->refoptionpckindex);
  gt_free(arguments);
}

static GtOptionParser *gt_idxlocali_option_parser_new(void *tool_arguments)
{
  IdxlocaliOptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionesaindex, *optionpckindex, *optiononline, *optioncmp;

  gt_assert(arguments != NULL);
  arguments->indexname = gt_str_new ();
  arguments->queryfiles = gt_str_array_new ();

  op = gt_option_parser_new
    ("[options] -q query-file-names [-esa|-pck] indexname",
     "Find all local alignments using suffix tree.");

  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");
  option = gt_option_new_filename_array("q","Specify files containing the "
                                            "query sequences",
                                        arguments->queryfiles);
  gt_option_parser_add_option (op, option);

  option = gt_option_new_long("match",
                              "Specify match score",
                              &arguments->matchscore, 1L);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_long("mismatch",
                              "Specify mismatch score",
                              &arguments->mismatchscore, -3L);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_long("gapstart",
                              "Specify gap start score",
                              &arguments->gapstart, -5L);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_long("gapextend",
                              "Specify gap extension score",
                              &arguments->gapextend, -2L);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong_min("th",
                                   "Specify the threshold",
                                    &arguments->threshold, 0, 1UL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  optionesaindex = gt_option_new_string("esa",
                                        "Specify index "
                                        "(enhanced suffix array)",
                                        arguments->indexname, NULL);
  gt_option_parser_add_option(op, optionesaindex);
  arguments->refoptionesaindex = gt_option_ref(optionesaindex);

  optionpckindex = gt_option_new_string("pck",
                                        "Specify index (packed index)",
                                        arguments->indexname, NULL);
  gt_option_parser_add_option(op, optionpckindex);
  arguments->refoptionpckindex = gt_option_ref (optionpckindex);
  gt_option_exclude (optionesaindex, optionpckindex);
  gt_option_is_mandatory_either(optionesaindex, optionpckindex);

  optiononline = gt_option_new_bool("online","Perform online searches",
                                    &arguments->doonline, false);
  gt_option_parser_add_option(op, optiononline);
  gt_option_is_development_option(optiononline);

  optioncmp = gt_option_new_bool("cmp","Compare results of offline and online "
                                 "searches",
                                 &arguments->docompare, false);
  gt_option_parser_add_option(op,optioncmp);
  gt_option_exclude(optiononline,optioncmp);

  option = gt_option_new_bool("s",
                              "Show alignments",
                              &arguments->showalignment, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_idxlocali_arguments_check(GT_UNUSED int rest_argc,
                                        void *tool_arguments,
                                        GT_UNUSED GtError * err)
{
  IdxlocaliOptions *arguments = tool_arguments;

  if (gt_option_is_set(arguments->refoptionesaindex))
  {
    arguments->withesa = true;
  } else
  {
    gt_assert(gt_option_is_set(arguments->refoptionpckindex));
    arguments->withesa = false;
  }
  return 0;
}

static int gt_idxlocali_runner(GT_UNUSED int argc,
                               GT_UNUSED const char **argv,
                               GT_UNUSED int parsed_args,
                               void *tool_arguments,
                               GtError * err)
{
  IdxlocaliOptions *arguments = tool_arguments;
  bool haserr = false;
  unsigned long idx;

  gt_error_check(err);
  gt_assert(arguments != NULL);

  gt_assert(parsed_args == argc);
  printf("# indexname(%s)=%s\n", arguments->withesa ? "esa" : "pck",
         gt_str_get(arguments->indexname));
  for (idx = 0; idx < gt_str_array_size (arguments->queryfiles); idx++)
  {
    printf("# queryfile=%s\n",gt_str_array_get (arguments->queryfiles, idx));
  }
  printf("# threshold=%lu\n", arguments->threshold);
  if (!haserr && gt_runidxlocali (arguments, err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}

GtTool *gt_idxlocali(void)
{
  return gt_tool_new(gt_idxlocali_arguments_new,
                     gt_idxlocali_arguments_delete,
                     gt_idxlocali_option_parser_new,
                     gt_idxlocali_arguments_check,
                     gt_idxlocali_runner);
}
