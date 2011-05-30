/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/cstr_array.h"
#include "core/error.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "extended/toolbox.h"
#include "match/cgr_spacedseed.h"
#include "tools/gt_congruence.h"

typedef struct
{
  bool withesa, docompare, verbose;
  GtStr *str_inputindex;
  GtStrArray *queryfilenames;
  GtOption *refoptionesaindex, *refoptionpckindex;
} Cge_spacedseed_options;

static void *gt_cge_spacedseed_arguments_new(void)
{
  Cge_spacedseed_options *arguments
    = gt_malloc(sizeof (Cge_spacedseed_options));
  arguments->str_inputindex = gt_str_new();
  arguments->queryfilenames = gt_str_array_new();
  return arguments;
}

static void gt_cge_spacedseed_arguments_delete(void *tool_arguments)
{
  Cge_spacedseed_options *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->str_inputindex);
  gt_str_array_delete(arguments->queryfilenames);
  gt_option_delete(arguments->refoptionpckindex);
  gt_option_delete(arguments->refoptionesaindex);
  gt_free(arguments);
}

static GtOptionParser
            *gt_cge_spacedseed_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option,
           *optionesaindex,
           *optionpckindex;
  Cge_spacedseed_options *arguments = tool_arguments;

  op = gt_option_parser_new("[options]",
                            "Match spaced seeds.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  optionesaindex = gt_option_new_string("esa",
                                     "Specify index (enhanced suffix array)",
                                     arguments->str_inputindex, NULL);
  gt_option_parser_add_option(op, optionesaindex);
  arguments->refoptionesaindex = gt_option_ref(optionesaindex);

  optionpckindex = gt_option_new_string("pck",
                                     "Specify index (packed index)",
                                     arguments->str_inputindex, NULL);
  gt_option_parser_add_option(op, optionpckindex);
  arguments->refoptionpckindex = gt_option_ref(optionpckindex);
  gt_option_exclude(optionesaindex,optionpckindex);
  gt_option_is_mandatory_either(optionesaindex,optionpckindex);

  option = gt_option_new_filename_array("q",
                                    "Specify files containing the "
                                    "query sequences",
                                    arguments->queryfilenames);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  option = gt_option_new_bool("cmp","compare results of offline and online "
                              "searches",&arguments->docompare,false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_cge_spacedseed_arguments_check(int rest_argc,
                                             void *tool_arguments,
                                             GtError *err)
{
  Cge_spacedseed_options *arguments = tool_arguments;

  if (gt_str_length(arguments->str_inputindex) == 0)
  {
    gt_error_set(err,"missing indexname");
    return -1;
  }
  if (gt_option_is_set(arguments->refoptionesaindex))
  {
    arguments->withesa = true;
  } else
  {
    gt_assert(gt_option_is_set(arguments->refoptionpckindex));
    arguments->withesa = false;
  }
  if (rest_argc != 0)
  {
    gt_error_set(err,"superfluous file arguments");
    return -1;
  }
  return 0;
}

static int gt_cge_spacedseed_runner(GT_UNUSED int argc,
                                    GT_UNUSED const char **argv,
                                    GT_UNUSED int parsed_args,
                                    void *tool_arguments,
                                    GtError *err)
{
  Cge_spacedseed_options *arguments = tool_arguments;
  GtLogger *logger = NULL;
  bool haserr = false;

  gt_assert(parsed_args == argc);
  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (arguments->verbose)
  {
    unsigned long idx;

    printf("# %sindex=%s\n",arguments->withesa ? "esa" : "pck",
                            gt_str_get(arguments->str_inputindex));
    for (idx = 0; idx < gt_str_array_size(arguments->queryfilenames); idx++)
    {
      printf("# queryfile=%s\n",
             gt_str_array_get(arguments->queryfilenames,idx));
    }
  }
  if (gt_matchspacedseed(arguments->withesa,
                      arguments->docompare,
                      gt_str_get(arguments->str_inputindex),
                      arguments->queryfilenames,
                      arguments->verbose,
                      err) != 0)
  {
    haserr = true;
  }
  gt_logger_delete(logger);
  return haserr ? - 1 : 0;
}

static GtTool* gt_cge_spacedseed(void)
{
  return gt_tool_new(gt_cge_spacedseed_arguments_new,
                     gt_cge_spacedseed_arguments_delete,
                     gt_cge_spacedseed_option_parser_new,
                     gt_cge_spacedseed_arguments_check,
                     gt_cge_spacedseed_runner);
}

static void *gt_cge_arguments_new(void)
{
  GtToolbox *cge_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(cge_toolbox, "spacedseed", gt_cge_spacedseed());
  return cge_toolbox;
}

static void gt_cge_arguments_delete(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  if (!index_toolbox) return;
  gt_toolbox_delete(index_toolbox);
}

static GtOptionParser* gt_cge_option_parser_new(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtOptionParser *op;

  gt_assert(index_toolbox != NULL);
  op = gt_option_parser_new(
                    "[option ...] congruence_tool [argument ...]",
                    "Call congruence tool with name congruence_tool and pass "
                    "argument(s) to it.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, index_toolbox);
  gt_option_parser_set_min_args(op, 1U);
  gt_option_parser_refer_to_manual(op);
  return op;
}

static int gt_cge_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  char **nargv = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(index_toolbox != NULL);

  /* determine tool */
  if (!gt_toolbox_has_tool(index_toolbox, argv[parsed_args]))
  {
    gt_error_set(err, "congruence tool '%s' not found; option -help lists "
                   "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call sub-tool */
  if (!had_err)
  {
    if (!(toolfunc = gt_toolbox_get(index_toolbox, argv[parsed_args])))
    {
      tool = gt_toolbox_get_tool(index_toolbox, argv[parsed_args]);
      gt_assert(tool != NULL);
    }
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    if (toolfunc != NULL)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv,
                            err);
  }
  gt_cstr_array_delete(nargv);
  return had_err;
}

GtTool* gt_congruence(void)
{
  return gt_tool_new(gt_cge_arguments_new,
                     gt_cge_arguments_delete,
                     gt_cge_option_parser_new,
                     NULL,
                     gt_cge_runner);
}
