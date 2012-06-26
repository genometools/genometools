/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#include <float.h>
#include <math.h>
#include <stdio.h>

#include "core/encseq_api.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/timer_api.h"
#include "core/showtime.h"
#include "core/str_array_api.h"
#include "core/unused_api.h"

#include "match/esa-map.h"
#include "match/esa-shulen.h"
#include "match/genomediff.h"
#include "match/idx-limdfs.h"
#include "match/sarr-def.h"
#include "match/shu-genomediff-pck-simple.h"
#include "match/shu-genomediff.h"

#include "tools/gt_genomediff.h"

static void* gt_genomediff_arguments_new(void)
{
  GtGenomediffArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->indexname = gt_str_new();
  arguments->unitfile = gt_str_new();
  arguments->queryname = gt_str_array_new();
  return arguments;
}

static void gt_genomediff_arguments_delete(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->indexname);
  gt_str_delete(arguments->unitfile);
  gt_str_array_delete(arguments->queryname);
  gt_option_delete(arguments->ref_esaindex);
  gt_option_delete(arguments->ref_pckindex);
  gt_option_delete(arguments->ref_queryname);
  gt_option_delete(arguments->ref_unitfile);
  gt_free(arguments);
}

static GtOptionParser* gt_genomediff_option_parser_new(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionquery, *optionesaindex, *optionpckindex,
           *optiontraverse, *optionscan, *option_unitfile;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -pck|-esa indexname ",
                            "Calculates shulens of genomes and the Kr.");

  /* -maxdepth */
  option =  gt_option_new_int("maxdepth", "max depth of .pbi-file",
                              &arguments->user_max_depth, -1);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -max_n */
  option = gt_option_new_ulong("max_n", "Number of precalculated values "
                             "for ln(n!) and pmax(x)",
                             &arguments->max_ln_n_fac, 1000UL);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* -esa */
  optionesaindex = gt_option_new_string("esa",
                                     "Specify index (enhanced suffix array)",
                                     arguments->indexname, NULL);
  gt_option_parser_add_option(op, optionesaindex);

  /* -pck */
  optionpckindex = gt_option_new_string("pck",
                                        "Specify index (packed index)",
                                        arguments->indexname, NULL);
  gt_option_parser_add_option(op, optionpckindex);

  gt_option_exclude(optionesaindex,optionpckindex);
  gt_option_is_mandatory_either(optionesaindex,optionpckindex);

  /* ref esa */
  arguments->ref_esaindex = gt_option_ref(optionesaindex);

  /* ref pck */
  arguments->ref_pckindex = gt_option_ref(optionpckindex);

  /*-unitfile*/
  option_unitfile =
    gt_option_new_filename("unitfile",
                           "specifies genomic units, "
                           "File is in lua-format like this:\n"
                           "units = {\n"
                           "  genome1 = { \"file1\", \"file2\" },\n"
                           "  genome2 = { \"file3\", \"file4\" }\n"
                           "}\n"
                           "only give basenames of the files! comment lines "
                           "start with '--' and will be ignored. See "
                           "GTDIR/testdata/genomediff/unitfile1.lua for an "
                           "example.",
                           arguments->unitfile);
  gt_option_parser_add_option(op, option_unitfile);

  /*ref unitfile*/
  arguments->ref_unitfile = gt_option_ref(option_unitfile);

  /* -query */
  optionquery = gt_option_new_filename_array("query",
                                       "Files containing the query sequences "
                                       "if this option is set a simple "
                                       "shustring search will be used." ,
                                       arguments->queryname);
  gt_option_is_extended_option(optionquery);
  gt_option_exclude(optionquery, option_unitfile);
  gt_option_parser_add_option(op, optionquery);

  /* ref query */
  arguments->ref_queryname = gt_option_ref(optionquery);

  /* thresholds */
  /* divergence error */
  option = gt_option_new_double("thr",
                                "Threshold for difference (du, dl) in "
                                "divergence calculation.\n"
                                "default: 1e-9",
                                &arguments->divergence_threshold,
                                pow(10.0, -9.0));
  gt_option_is_extended_option(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* expected shulen error */
  option = gt_option_new_double("abs_err",
                                "absolute error for expected shulen "
                                "calculation.\n"
                                "default: 1e-5",
                                &arguments->divergence_abs_err,
                                1e-5);
  gt_option_is_extended_option(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* relative expected shulen error */
  option = gt_option_new_double("rel_err",
                                "relative error for expected shulen "
                                "calculation.\n"
                                "default: 1e-3",
                                &arguments->divergence_rel_err,
                                1e-3);
  gt_option_is_extended_option(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* M */
  option = gt_option_new_double("M",
                                "threshold for minimum logarithm.\n"
                                "default: DBL_MIN",
                                &arguments->divergence_m,
                                DBL_MIN);
  gt_option_is_extended_option(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* shulen */
  option = gt_option_new_bool("shulen",
                              "prints sum of shulen and stops",
                              &arguments->shulen_only,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* scan */
  optionscan = gt_option_new_bool("scan",
                                  "do not load esa index but scan"
                                  " it sequentially",
                                  &arguments->scan,
                                  true);
  gt_option_exclude(optionscan, optionpckindex);
  gt_option_parser_add_option(op, optionscan);

  /* traverse */
  optiontraverse = gt_option_new_bool("traverse",
                              "traverses the virtual tree without calculating"
                              " anything, sets GT_ENV_OPTIONS=-showtime"
                              " does not work with -esa and -query",
                              &arguments->traverse_only,
                              false);
  gt_option_exclude(optiontraverse, optionesaindex);
  gt_option_exclude(optiontraverse, optionquery);
  gt_option_is_development_option(optiontraverse);
  gt_option_parser_add_option(op, optiontraverse);
  /* mail */
  gt_option_parser_set_mail_address(op, "<willrodt@zbh.uni-hamburg.de>");
  return op;
}

static int gt_genomediff_arguments_check(GT_UNUSED int rest_argc,
                                         void *tool_arguments,
                                         GT_UNUSED GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (gt_option_is_set(arguments->ref_esaindex))
  {
    arguments->with_esa = true;
  }
  else
  {
    gt_assert(gt_option_is_set(arguments->ref_pckindex));
    arguments->with_esa = false;
  }
  if (gt_option_is_set(arguments->ref_queryname))
  {
    arguments->simplesearch = true;
  }
  else
  {
    arguments->simplesearch = false;
  }
  if (!had_err && arguments->traverse_only &&
      !gt_showtime_enabled())
  {
    printf ("GT_ENV_OPTIONS should be set to -showtime\n");
    printf ("setting showtime = true\n");
    gt_showtime_enable();
  }
  if (!had_err && gt_option_is_set(arguments->ref_unitfile))
  {
    arguments->with_units = true;
  }
  else
  {
    arguments->with_units = false;
  }
  return had_err;
}

static int callpairswisesshulendist(const char *indexname,
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
    }
    else
    {
      printf("%lu\n",totalgmatchlength);
    }
  }
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static int gt_genomediff_runner(GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments, GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  int had_err = 0;
  GtLogger *logger;
  GtTimer *timer = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose,
                         GT_LOGGER_DEFLT_PREFIX,
                         stdout);
  gt_assert(logger);

  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description("start");
    gt_timer_start(timer);
    gt_assert(timer);
  }
  if (arguments->with_units)
  {
    gt_logger_log(logger, "unitfile option set, filename is %s\n",
       gt_str_get(arguments->unitfile));
  }
  if (!had_err)
  {
    if (arguments->simplesearch)
    {
      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "run simple search", stdout);
      }
      if (arguments->with_esa)
      {
        had_err = callpairswisesshulendist(gt_str_get(arguments->indexname),
                                           arguments->queryname,
                                           logger,
                                           err);
      }
      else
      {
        had_err = gt_genomediff_pck_shu_simple(logger,
                                               arguments,
                                               err);
      }
    }
    else
    {
      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "start shu search", stdout);
      }
      had_err = gt_genomediff_shu(logger,
                                  arguments,
                                  timer,
                                  err);
    }
  }
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  gt_logger_delete(logger);

  return had_err;
}

GtTool* gt_genomediff(void)
{
  return gt_tool_new(gt_genomediff_arguments_new,
                     gt_genomediff_arguments_delete,
                     gt_genomediff_option_parser_new,
                     gt_genomediff_arguments_check,
                     gt_genomediff_runner);
}
