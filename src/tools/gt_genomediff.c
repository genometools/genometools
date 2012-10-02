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
#include <string.h>

#include "core/array2dim_api.h"
#include "core/cstr_api.h"
#include "core/encseq.h"
#include "core/encseq_options.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/showtime.h"
#include "core/warning_api.h"
#include "extended/gtdatahelp.h"
#include "match/esa-fileend.h"
#include "match/genomediff_opt.h"
#include "match/index_options.h"
#include "match/sfx-opt.h"
#include "match/sfx-run.h"
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
  gt_encseq_options_delete(arguments->loadopts);
  gt_index_options_delete(arguments->idxopts);
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
  op = gt_option_parser_new("[option ...] "
                          "(INDEX | -indexname NAME SEQFILE SEQFILE [...]) ",
                          "Calculates Kr: pairwise distances between genomes.");

  /* options */
  option = gt_option_new_choice("indextype", "specify type of index, one of: "
                                "esa|pck|encseq. Where encseq is an encoded "
                                "sequence and an enhanced suffix array will be "
                                "constructed only in memory.",
                                arguments->indextype, indextypes[2],
                                indextypes);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("indexname", "Basename of encseq to construct.",
                                arguments->indexname, NULL);
  gt_option_parser_add_option(op, option);

  /*-unitfile*/
  option_unitfile =
    gt_option_new_filename("unitfile",
                           "specifies genomic units, "
                           "see below for description.",
                           arguments->unitfile);
  gt_option_parser_add_option(op, option_unitfile);
  arguments->ref_unitfile = gt_option_ref(option_unitfile);

  /* encseq options */
  arguments->loadopts =
    gt_encseq_options_register_loading(op, arguments->indexname);

  gt_option_is_development_option(
                        gt_encseq_options_lossless_option(arguments->loadopts));
  /* esa options */
  arguments->idxopts =
    gt_index_options_register_esa_noout(op);
  gt_option_is_development_option(
                            gt_index_options_spmopt_option(arguments->idxopts));

  /* scan */
  option = gt_option_new_bool("scan", "do not load esa index but scan "
                              "it sequentially.", &arguments->scan, true);
  gt_option_is_extended_option(option);
  gt_option_parser_add_option(op, option);

  /* dev options */
  /* -max_n */
  option = gt_option_new_ulong("max_n", "Number of precalculated values "
                               "for ln(n!) and pmax(x).",
                               &arguments->max_ln_n_fac, 1000UL);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -maxdepth */
  option =  gt_option_new_int("maxdepth", "max depth of .pbi-file, use with "
                              "-indextype pck.",
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

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* mail */
  gt_option_parser_set_mail_address(op, "<willrodt@zbh.uni-hamburg.de>");
  /* doc */
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
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

  if (rest_argc == 0) {
    gt_error_set(err, "give at least one file (base)name!");
    had_err = -1;
  }
  if (!had_err) {
    if (strcmp("esa", gt_str_get(arguments->indextype)) == 0)
      arguments->with_esa = true;
    else if (strcmp("pck", gt_str_get(arguments->indextype)) == 0)
      arguments->with_pck = true;
  }
  prepared_index = (arguments->with_esa || arguments->with_pck);

  if (!had_err && arguments->user_max_depth != -1 && !arguments->with_pck)
    gt_warning("option -maxdepth does only apply to -indextype pck");

  if (!had_err &&
      prepared_index && gt_encseq_options_mirrored_value(arguments->loadopts))
    gt_warning("option -mirrored is ignored with esa and pck index");

  if (!had_err && prepared_index && rest_argc > 1) {
    gt_error_set(err, "there should be only one basename argument with "
                 "-indextype esa|pck");
    had_err = -1;
  }
  if (rest_argc == 1 && gt_str_length(arguments->indexname) != 0) {
    gt_error_set(err, "Option -indexname is only needed with sequence files, "
                 "if one file is given as argument, this should be an index.");
    had_err = -1;
  }
  if (!had_err && rest_argc > 1 && gt_str_length(arguments->indexname) == 0) {
    gt_error_set(err, "use -indexname for basename of encseq");
    had_err = -1;
  }

  if (!had_err)
    arguments->with_units = gt_option_is_set(arguments->ref_unitfile);

  return had_err;
}

static int gt_genomediff_runner(int argc, const char **argv,
                                int parsed_args, void *tool_arguments,
                                GtError *err)
{
  bool mirrored = false;
  int had_err = 0,
      i;
  GtEncseq              *encseq = NULL;
  GtGenomediffArguments *arguments = tool_arguments;
  GtLogger              *logger;
  GtShuUnitFileInfo     *unit_info = NULL;
  GtTimer               *timer = NULL;

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

  if (gt_str_array_size(arguments->filenames) > 1UL) {
    GtEncseqEncoder *ee = gt_encseq_encoder_new();
    gt_encseq_encoder_set_timer(ee, timer);
    gt_encseq_encoder_set_logger(ee, logger);
    /* kr only makes sense for dna, so we can check this already with ee */
    gt_encseq_encoder_set_input_dna(ee);
    had_err = gt_encseq_encoder_encode(ee, arguments->filenames,
                                       gt_str_get(arguments->indexname), err);
    gt_encseq_encoder_delete(ee);
  }
  else {
    gt_str_append_str(arguments->indexname,
                      gt_str_array_get_str(arguments->filenames, 0));
    if (arguments->with_esa || arguments->with_pck) {
      GtStr *current_line = gt_str_new();
      FILE *prj_fp;
      const char *buffer;
      char **elements = NULL;

      prj_fp = gt_fa_fopen_with_suffix(gt_str_get(arguments->indexname),
                                       PROJECTFILESUFFIX,"rb",err);
      if (prj_fp == NULL)
        had_err = -1;
      while (!had_err && gt_str_read_next_line(current_line, prj_fp) != EOF) {
        buffer = gt_str_get(current_line);
        if (elements != NULL) {
          gt_free(elements[0]);
          gt_free(elements[1]);
        }
        gt_free(elements);
        elements = gt_cstr_split(buffer, '=');
        gt_log_log("%s", elements[0]);
        if (strcmp("mirrored", elements[0]) == 0) {
          gt_log_log("%s", elements[1]);
          if (strcmp("1", elements[1]) == 0) {
            mirrored = true;
            gt_log_log("sequences are treated as mirrored");
          }
        }
        gt_str_reset(current_line);
      }
      gt_str_delete(current_line);
      if (elements != NULL) {
        gt_free(elements[0]);
        gt_free(elements[1]);
      }
      gt_free(elements);
      gt_fa_xfclose(prj_fp);
    }
  }

  if (!had_err) {
    GtEncseqLoader *el = gt_encseq_loader_new_from_options(arguments->loadopts,
                                                           err);
    if (mirrored)
      gt_encseq_loader_mirror(el);
    encseq =
      gt_encseq_loader_load(el, gt_str_get(arguments->indexname), err);
    gt_encseq_loader_delete(el);
  }
  if (encseq == NULL)
    had_err = -1;
  if (!had_err) {
    unit_info = gt_shu_unit_info_new(encseq);
    if (arguments->with_units)
      had_err = gt_shu_unit_file_info_read(arguments->unitfile, unit_info,
                                           logger, err);
  }

  if (!had_err) {
    uint64_t **shusums = NULL;
    if (arguments->with_esa || arguments->with_pck) {
      shusums = gt_genomediff_shulen_sum(arguments, unit_info,
                                         logger, timer, err);
      if (shusums == NULL)
        had_err = -1;
    }
    else {
      const bool doesa = true;
      GenomediffInfo gd_info;
      Suffixeratoroptions sopts;
      sopts.beverbose = arguments->verbose;
      sopts.indexname = arguments->indexname;
      sopts.db = NULL;
      sopts.encopts = NULL;
      sopts.genomediff = true;
      sopts.inputindex = arguments->indexname;
      sopts.loadopts = arguments->loadopts;
      sopts.showprogress = false;
      sopts.idxopts = arguments->idxopts;

      gt_assert(unit_info != NULL);
      gt_array2dim_calloc(shusums, unit_info->num_of_genomes,
                          unit_info->num_of_genomes);
      gd_info.shulensums = shusums;
      gd_info.unit_info = unit_info;
      had_err = runsuffixerator(doesa, &sopts, &gd_info, logger, err);
    }
    if (!had_err && shusums != NULL) {
      had_err = gt_genomediff_kr_calc(shusums, arguments, unit_info,
                                      arguments->with_pck, logger, timer, err);
      gt_array2dim_delete(shusums);
    }
  }

  if (timer != NULL) {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  gt_logger_delete(logger);
  gt_encseq_delete(encseq);
  gt_shu_unit_info_delete(unit_info);

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
