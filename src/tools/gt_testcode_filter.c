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
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "extended/testcode_filter_visitor.h"
#include "extended/testcode_filter_stream.h"
#include "extended/visitor_stream.h"
#include "tools/gt_testcode_filter.h"

typedef struct {
  GtSeqid2FileInfo *s2fi;
  unsigned int windowsize;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  double threshold;
} GtTcodeFilterArguments;

static void* gt_testcode_filter_arguments_new(void)
{
  GtTcodeFilterArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_testcode_filter_arguments_delete(void *tool_arguments)
{
  GtTcodeFilterArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_testcode_filter_option_parser_new(void
                                                                *tool_arguments)
{
  GtTcodeFilterArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [GFF3_file]",
                            "Annotates CDS sequences with their mean "
                            "TESTCODE scores and filters out non-coding ones.");

  /* -seqfile, -matchdesc, -usedesc and -regionmapping */
  gt_seqid2file_register_options(op, arguments->s2fi);

  /* -width */
  option = gt_option_new_uint("windowsize",
                              "window size for testcode algorithm",
                              &arguments->windowsize, 24);
  gt_option_parser_add_option(op, option);

  /* -threshold */
  option = gt_option_new_double("threshold",
                                "threshold for coding classification",
                                &arguments->threshold, 0.95);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  return op;
}

static int gt_testcode_filter_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args, void *tool_arguments,
                                 GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL,
               *testcode_ann_stream = NULL,
               *testcode_filter_stream = NULL,
               *gff3_out_stream = NULL;
  GtTcodeFilterArguments *arguments = tool_arguments;
  GtRegionMapping *region_mapping = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (!had_err) {
    /* create gff3 input stream */
    gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);

    /* create region mapping */
    region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!region_mapping)
      had_err = -1;
  }

  if (!had_err) {
    testcode_ann_stream = gt_visitor_stream_new(gff3_in_stream,
                                gt_testcode_filter_visitor_new(region_mapping,
                                                        arguments->windowsize));
    testcode_filter_stream = gt_testcode_filter_stream_new(testcode_ann_stream,
                                                          arguments->threshold);
    gff3_out_stream = gt_gff3_out_stream_new(testcode_filter_stream,
                                             arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(gff3_out_stream, err);
  }

  /* free */
  gt_node_stream_delete(testcode_ann_stream);
  gt_node_stream_delete(testcode_filter_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(gff3_out_stream);
  gt_region_mapping_delete(region_mapping);

  return had_err;
}

GtTool* gt_testcode_filter(void)
{
  return gt_tool_new(gt_testcode_filter_arguments_new,
                     gt_testcode_filter_arguments_delete,
                     gt_testcode_filter_option_parser_new,
                     NULL,
                     gt_testcode_filter_runner);
}
