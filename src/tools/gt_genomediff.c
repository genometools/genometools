/*
  Copyright (c) 2010 Willrodt <dwillrodt@zbh.uni-hamburg.de>
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
#include "core/str_array_api.h"
#include "core/unused_api.h"

#include "match/genomediff.h"
#include "match/idx-limdfs.h"
#include "match/shu-genomediff-kr2.h"
#include "match/shu-genomediff-simple.h"

#include "tools/gt_genomediff.h"

static void* gt_genomediff_arguments_new(void)
{
  GtGenomediffArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->indexname = gt_str_new();
  arguments->queryname = gt_str_array_new();
  return arguments;
}

static void gt_genomediff_arguments_delete(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->indexname);
  gt_str_array_delete(arguments->queryname);
  gt_option_delete(arguments->ref_esaindex);
  gt_option_delete(arguments->ref_pckindex);
  gt_option_delete(arguments->ref_queryname);
  gt_free(arguments);
}

static GtOptionParser* gt_genomediff_option_parser_new(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionquery, *optionesaindex, *optionpckindex;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -pck indexname "
                            "[-query sequencefile(s)]",
                            "Reads in a packedindex.");

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
  gt_option_is_development_option(optionesaindex);
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

  /* -query */
  optionquery = gt_option_new_filenamearray("query",
                                       "Files containing the query sequences "
                                       "if this option is set a simple "
                                       "shustring search will be used." ,
                                       arguments->queryname);
  gt_option_parser_add_option(op, optionquery);

  /* ref query */
  arguments->ref_queryname = gt_option_ref(optionquery);

  /* thresholds */
  option = gt_option_new_double("thr",
                                "Threshold for difference (du, dl) in"
                                "divergence calculation.\n"
                                "default: 1e-9",
                                &arguments->divergence_threshold,
                                pow(10.0, -9.0));
  gt_option_is_development_option(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_double("abs_err",
                                "absolut error for epected shulen "
                                "calculation.\n"
                                "default: 1e-5",
                                &arguments->divergence_abs_err,
                                1e-5);
  gt_option_is_development_option(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_double("rel_err",
                                "relative error for epected shulen "
                                "calculation.\n"
                                "default: 1e-3",
                                &arguments->divergence_rel_err,
                                1e-3);
  gt_option_is_development_option(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_double("M",
                                "threshold for minimum logarithm.\n"
                                "default: DBL_MIN",
                                &arguments->divergence_m,
                                DBL_MIN);
  gt_option_is_development_option(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("shulen",
                              "prints sum of shulen and stops",
                              &arguments->shulen_only,
                              false);
  gt_option_parser_add_option(op, option);

  /* mail */
  gt_option_parser_set_mailaddress(op, "<dwillrodt@zbh.uni-hamburg.de>");
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

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). */
  if (gt_option_is_set(arguments->ref_esaindex))
  {
    arguments->withesa = true;
  } else
  {
    gt_assert(gt_option_is_set(arguments->ref_pckindex));
    arguments->withesa = false;
  }
  if (gt_option_is_set(arguments->ref_queryname))
    arguments->simplesearch = true;
  else
    arguments->simplesearch = false;
  if (arguments->withesa)
  {
    printf("not implemented option -esa used, sorry, try -pck instead\n");
    had_err = 1;
  }
  return had_err;
}

static int gt_genomediff_runner(GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments, GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  int had_err = 0;
  Genericindex *genericindexSubject;
  GtLogger *logger;
  const GtEncseq *encseq = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose,
                         GT_LOGGER_DEFLT_PREFIX,
                         stdout);
  gt_assert(logger);

  genericindexSubject = genericindex_new(gt_str_get(
                                           arguments->indexname),
                                         arguments->withesa,
                                         true,
                                         false,
                                         true,
                                         arguments->user_max_depth,
                                         logger,
                                         err);
  if (genericindexSubject == NULL)
    had_err = 1;
  else
    encseq = genericindex_getencseq(genericindexSubject);

  if (!had_err)
  {
    if (arguments->simplesearch)
      had_err = gt_genomediff_run_simple_search(genericindexSubject,
                                                encseq,
                                                logger,
                                                arguments,
                                                err);
    else
      had_err = gt_genomediff_run_kr2_search(genericindexSubject,
                                             encseq,
                                             logger,
                                             arguments,
                                             err);
  }

/*XXX*/

  genericindex_delete(genericindexSubject);
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
