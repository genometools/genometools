/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "core/tokenizer.h"
#include "core/io.h"
#include "extended/gtdatahelp.h"
#include "match/esa-fileend.h"
#include "match/genomediff_opt.h"
#include "match/index_options.h"
#include "match/sfx-opt.h"
#include "match/sfx-run.h"
#include "match/shu-genomediff.h"
#include "tools/gt_gdiffcalc.h"

static void* gt_gdiffcalc_arguments_new(void)
{
  GtGenomediffArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->indexname = gt_str_new();
  arguments->unitfile = gt_str_new();
  arguments->indextype = gt_str_new();
  arguments->filenames = gt_str_array_new();
  arguments->with_esa = arguments->with_pck = false;
  return arguments;
}

static void gt_gdiffcalc_arguments_delete(void *tool_arguments)
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

static GtOptionParser* gt_gdiffcalc_option_parser_new(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *option_unitfile;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] "
                          "-indexname NAME AVGSHULEN) ",
                          "Calculates Kr: pairwise distances between genomes.");

  /* options */
  option = gt_option_new_string("indexname", "Basename of encseq to construct.",
                                arguments->indexname, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /*-unitfile*/
  option_unitfile =
    gt_option_new_filename("unitfile",
                           "specifies genomic units, see below for description",
                           arguments->unitfile);
  gt_option_parser_add_option(op, option_unitfile);
  arguments->ref_unitfile = gt_option_ref(option_unitfile);

  /* encseq options */
  arguments->loadopts =
    gt_encseq_options_register_loading(op, arguments->indexname);

  gt_option_is_development_option(
                        gt_encseq_options_lossless_option(arguments->loadopts));

  /* dev options */
  /* -max_n */
  option = gt_option_new_ulong("max_n", "Number of precalculated values "
                               "for ln(n!) and pmax(x)",
                               &arguments->max_ln_n_fac, 1000UL);
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

static int gt_gdiffcalc_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  gt_assert(rest_argc == 1);
  if (!had_err)
    arguments->with_units = gt_option_is_set(arguments->ref_unitfile);

  return had_err;
}

static int gt_gdiffcalc_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  int had_err = 0, i;
  unsigned long lcounter = 0, zcounter = 0;
  double **shusums = NULL;
  GtEncseq              *encseq = NULL;
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
    timer = gt_timer_new_with_progress_description("load encseq");
    gt_timer_start(timer);
    gt_assert(timer);
  }

  if (arguments->with_units) {
    gt_logger_log(logger, "unitfile option set, filename is %s\n",
                  gt_str_get(arguments->unitfile));
  }

  if (!had_err) {
    GtEncseqLoader *el = gt_encseq_loader_new_from_options(arguments->loadopts,
                                                           err);
    encseq =
      gt_encseq_loader_load(el, gt_str_get(arguments->indexname), err);
    gt_encseq_loader_delete(el);
  }
  if (encseq == NULL)
    had_err = -1;

  if (timer != NULL)
    gt_timer_show_progress(timer, "load units", stdout);

  if (!had_err) {
    unit_info = gt_shu_unit_info_new(encseq);
    if (arguments->with_units)
      had_err = gt_shu_unit_file_info_read(arguments->unitfile, unit_info,
                                           logger, err);
  }

  if (timer != NULL)
    gt_timer_show_progress(timer, "read table", stdout);

  if (!had_err) {
    GtIO *table_file = NULL;
    GtTokenizer *tokenizer = NULL;
    GtStr *line = NULL;

    gt_assert(unit_info != NULL);
    gt_array2dim_calloc(shusums, unit_info->num_of_genomes,
                        unit_info->num_of_genomes);

    table_file = gt_io_new(gt_str_array_get(arguments->filenames, 0), "r");
    tokenizer = gt_tokenizer_new(table_file);
    line = gt_tokenizer_get_token(tokenizer);
    while (line != NULL && !had_err) {
      char *cline = gt_str_get(line);
      char *elem = strtok(cline, ";");
      zcounter = 0;
      while (elem != NULL && !had_err) {
        if (*elem != '#') {
          if (1 != sscanf(elem, "%lf",
                          &shusums[lcounter][zcounter])) {
            had_err = 1;
            gt_error_set(err, "couldn't scan");
            break;
          }
          gt_logger_log(logger,"wert: %lf", shusums[lcounter][zcounter]);
          zcounter++;
        }
        else {
          gt_logger_log(logger, "name: %s", elem++);
        }
        elem = strtok(NULL, ";");
      }
      gt_tokenizer_next_token(tokenizer);
      gt_str_delete(line);
      line = gt_tokenizer_get_token(tokenizer);
      lcounter++;
      gt_logger_log(logger, "line %ld", lcounter);
    }
  }
  if (!had_err) {
    unsigned long num_of_seq, file_idx, seq_idx, startpos;
    GT_UNUSED unsigned long oldpos = 0;

    gt_assert(unit_info != NULL);
    gt_assert(lcounter == zcounter);
    gt_assert(lcounter == unit_info->num_of_genomes);

    num_of_seq = gt_encseq_num_of_sequences(unit_info->encseq);

    for (seq_idx = 0; seq_idx < num_of_seq; seq_idx++) {
      startpos = gt_encseq_seqstartpos(unit_info->encseq, seq_idx);
      file_idx = gt_encseq_filenum(unit_info->encseq, startpos);
      gt_log_log("seq: %lu starts at: %lu\n"
                 "belonges to file: %lu which is part of genome: %s",
                 seq_idx, startpos, file_idx,
                 gt_str_array_get(unit_info->genome_names,
                                  unit_info->map_files[file_idx]));
      gt_assert(oldpos <= startpos);
      oldpos = startpos;
    }
  }
  if (!had_err && shusums != NULL) {
    had_err = gt_genomediff_calculate_div_from_avg(shusums, arguments,
                                                   unit_info,
                                                   logger, timer, err);
    gt_array2dim_delete(shusums);
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

GtTool* gt_gdiffcalc(void)
{
  return gt_tool_new(gt_gdiffcalc_arguments_new,
                  gt_gdiffcalc_arguments_delete,
                  gt_gdiffcalc_option_parser_new,
                  gt_gdiffcalc_arguments_check,
                  gt_gdiffcalc_runner);
}
