/*
  Copyright (c) 2007-2011, 2013 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008       Center for Bioinformatics, University of Hamburg

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

#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "extended/chseqids_stream.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtdatahelp.h"
#include "extended/sort_stream_api.h"
#include "tools/gt_chseqids.h"

#define DEFAULT_JOINLENGTH 300

typedef struct {
  bool sort,
       verbose;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} ChseqidsArguments;

static void* gt_chseqids_arguments_new(void)
{
  ChseqidsArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_chseqids_arguments_delete(void *tool_arguments)
{
  ChseqidsArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_chseqids_option_parser_new(void *tool_arguments)
{
  ChseqidsArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] mapping_file [GFF3_file]",
                         "Change sequence ids by the mapping given in a "
                         "mapping file.");

  /* -sort */
  option = gt_option_new_bool("sort",
                              "sort the GFF3 features after changing the "
                              "sequence ids\n(memory consumption is "
                              "proportional to the input file size)",
                              &arguments->sort, false);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  /* parse options */
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_min_max_args(op, 1, 2);

  return op;
}

static int gt_chseqids_runnter(GT_UNUSED int argc, const char **argv,
                               int parsed_args, void *tool_arguments,
                               GtError *err)
{
  ChseqidsArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream, *chseqids_stream, *sort_stream = NULL,
               *gff3_out_stream = NULL;
  GtStr *chseqids;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* create the streams */
  gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args + 1]);
  if (arguments->verbose && arguments->outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);
  chseqids = gt_str_new_cstr(argv[parsed_args]);
  chseqids_stream = gt_chseqids_stream_new(gff3_in_stream, chseqids, err);
  if (!chseqids_stream)
    had_err = -1;
  gt_str_delete(chseqids);
  if (!had_err) {
    if (arguments->sort) {
      sort_stream = gt_sort_stream_new(chseqids_stream);
      gff3_out_stream = gt_gff3_out_stream_new(sort_stream, arguments->outfp);
    }
    else {
      gff3_out_stream = gt_gff3_out_stream_new(chseqids_stream,
                                               arguments->outfp);
    }
  }

  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(chseqids_stream);
  gt_node_stream_delete(sort_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_chseqids(void)
{
  return gt_tool_new(gt_chseqids_arguments_new,
                     gt_chseqids_arguments_delete,
                     gt_chseqids_option_parser_new,
                     NULL,
                     gt_chseqids_runnter);
}
