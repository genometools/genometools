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
#include "core/showtime.h"
#include "core/str_array_api.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "match/esa-map.h"
#include "match/esa-shulen.h"
#include "match/genomediff_opt.h"
#include "match/idx-limdfs.h"
#include "match/sarr-def.h"
#include "match/shu-genomediff.h"
#include "tools/gt_genomediff.h"

static void* gt_genomediff_arguments_new(void)
{
  GtGenomediffArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->indexname = gt_str_new();
  arguments->unitfile = gt_str_new();
  arguments->indextype = gt_str_new();
  arguments->filenames = gt_str_array_new();
  arguments->with_esa = arguments->with_pck = false;
  return arguments;
}

static void gt_genomediff_arguments_delete(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->indexname);
  gt_str_delete(arguments->unitfile);
  gt_str_delete(arguments->indextype);
  gt_str_array_delete(arguments->filenames);
  gt_option_delete(arguments->ref_unitfile);
  gt_free(arguments);
}

static GtOptionParser* gt_genomediff_option_parser_new(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *option_unitfile;
  static const char *indextypes[] = { "esa", "pck", "encseq", NULL };

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] (indexname|filenames) ",
                          "Calculates Kr: pairwise distances between genomes.");

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* -index */
  option = gt_option_new_choice("indextype", "specify type of index, one of: "
                                "esa|pck|encseq. Where enc is an encoded "
                                "sequence and an enhanced suffix array will be "
                                "constructed temporarily. Will be ignored if "
                                "sequence files are given.",
                                arguments->indextype, indextypes[2],
                                indextypes);

  gt_option_parser_add_option(op, option);

  /* -indexname */
  option = gt_option_new_string("indexname", "if files given are sequence "
                                "files, this will be the basename of the "
                                "encoded sequence created during index "
                                "construction.", arguments->indexname, NULL);
  gt_option_parser_add_option(op, option);

  /* -mirror */
  option = gt_option_new_bool("mirror", "specify if encseq should be mirrored "
                              "to include the reverse complement. Ignored with "
                              "esa and pck index",
                              &arguments->mirror, false);
  gt_option_parser_add_option(op, option);

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
  arguments->ref_unitfile = gt_option_ref(option_unitfile);

  /* scan */
  option = gt_option_new_bool("scan", "do not load esa index but scan "
                              "it sequentially", &arguments->scan, true);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* dev options */
  /* -max_n */
  option = gt_option_new_ulong("max_n", "Number of precalculated values "
                             "for ln(n!) and pmax(x)",
                             &arguments->max_ln_n_fac, 1000UL);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -maxdepth */
  option =  gt_option_new_int("maxdepth", "max depth of .pbi-file, use with "
                              "-indextype pck",
                              &arguments->user_max_depth, -1);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* thresholds */
  /* divergence error */
  option = gt_option_new_double("thr",
                                "Threshold for difference (du, dl) in "
                                "divergence calculation.\n"
                                "default: 1e-9",
                                &arguments->divergence_threshold,
                                1e-9);
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

  /* mail */
  gt_option_parser_set_mail_address(op, "<willrodt@zbh.uni-hamburg.de>");
  return op;
}

static int gt_genomediff_arguments_check(int rest_argc,
                                         void *tool_arguments,
                                         GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  bool prepared_index;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (!had_err) {
    if (strcmp("esa", gt_str_get(arguments->indextype)) == 0)
      arguments->with_esa = true;
    else if (strcmp("pck", gt_str_get(arguments->indextype)) == 0)
      arguments->with_pck = true;
  }
  prepared_index = (arguments->with_esa || arguments->with_pck);

  if (arguments->user_max_depth != -1 && !arguments->with_pck)
    gt_warning("option -maxdepth does only apply to -indextype pck");

  if (prepared_index && arguments->mirror)
    gt_warning("option -mirror is ignored with esa and pck index");

  if (rest_argc == 0) {
    gt_error_set(err, "at least one filename has to be given");
    had_err = -1;
  }
  if (rest_argc > 1) {
    if (gt_str_length(arguments->indexname) == 0) {
      had_err = -1;
      gt_error_set(err, "more than one file given, pleas give a basename for "
                   "the encseq to be created");
    }
    else if (prepared_index) {
      gt_error_set(err, "indextypes esa and pck only accept one filename "
                        "(basename)!");
      had_err = -1;
    }
  }

  if (!had_err)
    arguments->with_units = gt_option_is_set(arguments->ref_unitfile);

  return had_err;
}

static int gt_genomediff_runner(GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments, GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  int had_err = 0,
      i;
  GtLogger *logger;
  GtTimer *timer = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose,
                         GT_LOGGER_DEFLT_PREFIX,
                         stdout);
  gt_assert(logger);

  for (i = parsed_args; i < argc; i++) {
    gt_str_array_add_cstr(arguments->filenames, argv[i]);
  }

  if (gt_showtime_enabled()) {
    timer = gt_timer_new_with_progress_description("start");
    gt_timer_start(timer);
    gt_assert(timer);
  }

  if (arguments->with_units) {
    gt_logger_log(logger, "unitfile option set, filename is %s\n",
                  gt_str_get(arguments->unitfile));
  }

  if (timer != NULL)
    gt_timer_show_progress(timer, "start shu search", stdout);

  if (arguments->with_esa || arguments->with_pck)
    had_err = gt_genomediff_shulen_sum(logger, arguments, timer, err);
  else {
    had_err = -1;
    gt_error_set(err, "functionality not implemented yet");
  }

  if (timer != NULL) {
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
