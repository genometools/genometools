/*
  Copyright (c) 2007-2009, 2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "extended/genome_node.h"
#include "extended/gff3_in_stream_api.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtdatahelp.h"
#include "extended/uniq_stream_api.h"
#include "tools/gt_uniq.h"

typedef struct {
  bool verbose;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} UniqArguments;

static void* gt_uniq_arguments_new(void)
{
  UniqArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_uniq_arguments_delete(void *tool_arguments)
{
  UniqArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_uniq_option_parser_new(void *tool_arguments)
{
  UniqArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file]", "Filter out repeated "
                            "feature node graphs in a sorted GFF3 file.");

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_max_args(op, 1);

  return op;
}

static int gt_uniq_runner(GT_UNUSED int argc, const char **argv,
                          int parsed_args, void *tool_arguments, GtError *err)
{
  UniqArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream,
               *uniq_stream = NULL,
               *gff3_out_stream = NULL;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  /* create gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);
  if (arguments->verbose && arguments->outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create uniq stream */
  uniq_stream = gt_uniq_stream_new(gff3_in_stream);

  /* create gff3 output stream */
  gff3_out_stream = gt_gff3_out_stream_new(uniq_stream, arguments->outfp);

  /* pull the features through the stream and free them afterwards */
  had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(uniq_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_uniq(void)
{
  return gt_tool_new(gt_uniq_arguments_new,
                     gt_uniq_arguments_delete,
                     gt_uniq_option_parser_new,
                     NULL,
                     gt_uniq_runner);
}
