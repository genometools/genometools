/*
  Copyright (c) 2014-2016 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014-2016 Genome Research Ltd.

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

#include <unistd.h>
#include "core/cstr_api.h"
#include "core/error.h"
#include "core/fileutils_api.h"
#include "core/gtdatapath.h"
#include "core/ma.h"
#include "core/output_file_api.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/array_in_stream_api.h"
#include "extended/array_out_stream_api.h"
#include "extended/feature_stream_api.h"
#include "extended/feature_index_api.h"
#include "extended/feature_index_memory.h"
#include "extended/gff3_in_stream.h"
#include "extended/seqid2file.h"
#include "extended/sort_stream_api.h"
#include "extended/spec_visitor.h"
#include "extended/typecheck_info.h"
#include "extended/visitor_stream.h"
#include "tools/gt_speck.h"

typedef struct {
  GtStr *specfile, *format;
  bool verbose,
       colored,
       fail_hard,
       provideindex,
       sort;
  GtSeqid2FileInfo *s2fi;
  GtOutputFileInfo *ofi;
  GtTypecheckInfo *tci;
  GtFile *outfp;
} SpeccheckArguments;

static void *gt_speck_arguments_new(void)
{
  SpeccheckArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->specfile = gt_str_new();
  arguments->format = gt_str_new();
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  arguments->tci = gt_typecheck_info_new();
  arguments->outfp = NULL;
  return arguments;
}

static void gt_speck_arguments_delete(void *tool_arguments)
{
  SpeccheckArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->specfile);
  gt_str_delete(arguments->format);
  gt_file_delete(arguments->outfp);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_output_file_info_delete(arguments->ofi);
  gt_typecheck_info_delete(arguments->tci);
  gt_free(arguments);
}

static GtOptionParser* gt_speck_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option;
  SpeccheckArguments *arguments = tool_arguments;

  /* init */
  op = gt_option_parser_new("[options] [GFF3_file ...]",
                            "Checks spec definition compliance in GFF3 input.");

  option = gt_option_new_filename("specfile",
                                  "file with specification definition",
                                  arguments->specfile);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  option = gt_option_new_bool("colored", "show colored output",
                              &arguments->colored, true);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("provideindex", "provide feature index in "
                              "specfile namespace (requires O(n) memory for n "
                              "input features)",
                              &arguments->provideindex, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("sort", "sort input before checking (requires "
                              "O(n) memory for n input features)",
                              &arguments->sort, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("failhard", "stop processing and report runtime "
                              "errors instead of recording them in the results",
                              &arguments->fail_hard, false);
  gt_option_parser_add_option(op, option);

  /* -format */
  option = gt_option_new_string("output", "output format\n"
                                "choose from: [json, text, html, statsonly, "
                                "tabular] or give path to output driver",
                                arguments->format, "text");
  gt_option_parser_add_option(op, option);

  gt_typecheck_info_register_options_with_default(arguments->tci, op, "so");
  gt_seqid2file_register_options_ext(op, arguments->s2fi, false, false);
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_speck_arguments_check(GT_UNUSED int rest_argc,
                                    void *tool_arguments,
                                    GT_UNUSED GtError *err)
{
  SpeccheckArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if ((arguments->outfp || (!arguments->outfp && !isatty(STDOUT_FILENO)))
      && arguments->colored) {
    gt_warning("not printing to terminal, disabling colored output");
    arguments->colored = false;
  }

  return had_err = 0;
}

static void gt_speck_record_warning(void *data, const char *format, va_list ap)
{
  char buf[BUFSIZ];
  gt_assert(data && format);
  (void) vsnprintf(buf, (BUFSIZ * sizeof (char)), format, ap);
  gt_spec_results_record_warning((GtSpecResults*) data, buf);
}

static int gt_speck_runner(int argc, const char **argv, int parsed_args,
                               void *tool_arguments, GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL, *checker_stream = NULL,
               *a_in_stream = NULL, *a_out_stream = NULL,
               *feature_stream = NULL, *sort_stream = NULL,
               *last_stream = NULL;
  GtNodeVisitor *spec_visitor = NULL;
  GtSpecResults *res = NULL;
  GtFeatureIndex *fi = NULL;
  GtTypeChecker *type_checker = NULL;
  GtTimer *t = NULL;
  GtRegionMapping *rm = NULL;
  GtArray *arr = gt_array_new(sizeof (GtFeatureNode*));
  GtStr *prog, *speclib;
  SpeccheckArguments *arguments = tool_arguments;

  int had_err = 0;
  gt_error_check(err);

  res = gt_spec_results_new();
  gt_assert(res);

  if (gt_file_exists(gt_str_get(arguments->format))) {
    speclib = gt_str_ref(arguments->format);
  } else {
    prog = gt_str_new();
    gt_str_append_cstr_nt(prog, gt_error_get_progname(err),
                    gt_cstr_length_up_to_char(gt_error_get_progname(err), ' '));
    speclib = gt_get_gtdata_path(gt_str_get(prog), NULL);
    gt_str_delete(prog);
    gt_str_append_cstr(speclib, "/spec/output_drivers/");
    gt_str_append_str(speclib, arguments->format);

    if (!gt_file_exists(gt_str_get(speclib))) {
      gt_error_set(err, "output driver file \"%s\" does not exist",
                   gt_str_get(speclib));
      had_err = -1;
    }
  }

  if (!had_err) {
    spec_visitor = gt_spec_visitor_new(gt_str_get(arguments->specfile), res,
                                       err);
    if (!spec_visitor) {
      gt_spec_results_delete(res);
      return -1;
    }
  }

  t = gt_timer_new();
  gt_assert(t);

  /* add region mapping if given */
  if (!had_err && gt_seqid2file_option_used(arguments->s2fi)) {
    rm = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!rm)
      had_err = -1;
    if (!had_err)
      gt_spec_visitor_add_region_mapping((GtSpecVisitor*) spec_visitor, rm);
  }

  /* set type checker if necessary */
  if (!had_err && gt_typecheck_info_option_used(arguments->tci)) {
    type_checker = gt_typecheck_info_create_type_checker(arguments->tci, err);
    if (!type_checker)
      had_err = -1;
    if (!had_err)
      gt_spec_visitor_add_type_checker((GtSpecVisitor*) spec_visitor,
                                       type_checker);
  }

  if (!had_err) {
    /* set runtime error behaviour */
    if (arguments->fail_hard)
      gt_spec_visitor_fail_on_runtime_error((GtSpecVisitor*) spec_visitor);
    else
      gt_spec_visitor_report_runtime_errors((GtSpecVisitor*) spec_visitor);

    /* redirect warnings */
    gt_warning_set_handler(gt_speck_record_warning, res);

    last_stream = gff3_in_stream = gt_gff3_in_stream_new_unsorted(
                                                            argc - parsed_args,
                                                            argv + parsed_args);
    gt_assert(gff3_in_stream);
    gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream*) gff3_in_stream);

    /* insert sort stream if requested */
    if (arguments->sort) {
      last_stream = sort_stream = gt_sort_stream_new(last_stream);
    }

    /* if -provideindex is given, collect input features and index them first */
    if (arguments->provideindex) {
      fi = gt_feature_index_memory_new();
      gt_assert(fi);

      last_stream = feature_stream = gt_feature_stream_new(last_stream, fi);
      gt_assert(feature_stream);

      last_stream = a_out_stream = gt_array_out_stream_all_new(last_stream, arr,
                                                               err);
      if (!a_out_stream)
        had_err = -1;

      gt_timer_start(t);

      if (!had_err)
        had_err = gt_node_stream_pull(last_stream, err);

      if (!had_err) {
        gt_spec_visitor_add_feature_index((GtSpecVisitor*) spec_visitor,
                                          gt_feature_index_ref(fi));
        last_stream = a_in_stream = gt_array_in_stream_new(arr, NULL, err);
        if (!a_in_stream)
          had_err = -1;
      }
    } else {
      gt_timer_start(t);
    }

    if (!had_err) {
      checker_stream = gt_visitor_stream_new(last_stream, spec_visitor);
      gt_assert(checker_stream);
    }

    /* perform checking  */
    if (!had_err)
      had_err = gt_node_stream_pull(checker_stream, err);

    gt_timer_stop(t);

    /* reset warnings output */
    gt_warning_set_handler(gt_warning_default_handler, NULL);

    /* output results */
    if (!had_err) {
      GtStr *runtime = gt_str_new();
      gt_timer_get_formatted(t, GT_WD ".%06ld", runtime);
      had_err = gt_spec_results_render_template(res, gt_str_get(speclib),
                                                arguments->outfp,
                                                gt_str_get(arguments->specfile),
                                                arguments->verbose,
                                                arguments->colored,
                                                gt_str_get(runtime), err);
      gt_str_delete(runtime);
    }

    if (!had_err && (gt_spec_results_has_runtime_errors(res)
          || gt_spec_results_has_failures(res))) {
      had_err = -2;
    }
  }

  /* free */
  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(a_in_stream);
  gt_node_stream_delete(a_out_stream);
  gt_node_stream_delete(checker_stream);
  gt_node_stream_delete(feature_stream);
  gt_node_stream_delete(sort_stream);
  gt_spec_results_delete(res);
  gt_feature_index_delete(fi);
  gt_type_checker_delete(type_checker);
  gt_timer_delete(t);
  gt_array_delete(arr);
  gt_str_delete(speclib);

  return had_err;
}

GtTool* gt_speck(void)
{
  return gt_tool_new(gt_speck_arguments_new,
                     gt_speck_arguments_delete,
                     gt_speck_option_parser_new,
                     gt_speck_arguments_check,
                     gt_speck_runner);
}
