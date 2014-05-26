/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd.

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

#include "core/ma.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "extended/spec_visitor.h"
#include "extended/gff3_in_stream.h"
#include "extended/visitor_stream.h"
#include "tools/gt_speck.h"

typedef struct {
  GtStr *specfile;
  bool verbose,
       colored,
       fail_hard;
} SpeccheckArguments;

static void *gt_speck_arguments_new(void)
{
  SpeccheckArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->specfile = gt_str_new();
  return arguments;
}

static void gt_speck_arguments_delete(void *tool_arguments)
{
  SpeccheckArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->specfile);
  gt_free(arguments);
}

static GtOptionParser* gt_speck_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option;
  SpeccheckArguments *arguments = tool_arguments;

  /* init */
  op = gt_option_parser_new("[GFF3_file ...]",
                            "Checks spec definition compliance in GFF3 input.");

  option = gt_option_new_filename("specfile",
                                  "file with specification definition",
                                  arguments->specfile);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  option = gt_option_new_bool("colored", "show colored output",
                              &arguments->colored, true);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("failhard", "stop processing and report runtime "
                              "errors instead of recording them in the results",
                              &arguments->fail_hard, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_speck_runner(int argc, const char **argv, int parsed_args,
                               void *tool_arguments, GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL, *checker_stream = NULL;
  GtNodeVisitor *spec_visitor = NULL;
  GtSpecResults *res = NULL;
  GtTimer *t = NULL;
  SpeccheckArguments *arguments = tool_arguments;

  int had_err = 0;
  gt_error_check(err);

  res = gt_spec_results_new();
  spec_visitor = gt_spec_visitor_new(gt_str_get(arguments->specfile), res,
                                     err);
  if (!spec_visitor)
    return -1;

  if (arguments->fail_hard)
    gt_spec_visitor_fail_on_runtime_error((GtSpecVisitor*) spec_visitor);
  else
    gt_spec_visitor_report_runtime_errors((GtSpecVisitor*) spec_visitor);

  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  gt_assert(gff3_in_stream);

  checker_stream = gt_visitor_stream_new(gff3_in_stream, spec_visitor);
  gt_assert(checker_stream);

  t = gt_timer_new();
  gt_timer_start(t);
  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(checker_stream, err);
  gt_timer_stop(t);

  if (!had_err)
    gt_spec_results_report(res, NULL, gt_str_get(arguments->specfile),
                           arguments->verbose, arguments->colored);

  if (!had_err)
    gt_timer_show_formatted(t, "Finished in " GT_WD ".%06ld s.\n", stderr);

  /* free */
  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(checker_stream);
  gt_spec_results_delete(res);
  gt_timer_delete(t);

  return had_err;
}

GtTool* gt_speck(void)
{
  return gt_tool_new(gt_speck_arguments_new,
                     gt_speck_arguments_delete,
                     gt_speck_option_parser_new,
                     NULL,
                     gt_speck_runner);
}
