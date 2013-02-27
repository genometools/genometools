/*
  Copyright (c) 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/output_file_api.h"
#include "extended/feature_type.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/inter_feature_stream_api.h"
#include "tools/gt_interfeat.h"

typedef struct {
  GtStr *outside_type,
        *inter_type;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} InterFeatArguments;

static void* gt_interfeat_arguments_new(void)
{
  InterFeatArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->outside_type = gt_str_new();
  arguments->inter_type = gt_str_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_interfeat_arguments_delete(void *tool_arguments)
{
  InterFeatArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_delete(arguments->inter_type);
  gt_str_delete(arguments->outside_type);
  gt_free(arguments);
}

static GtOptionParser* gt_interfeat_option_parser_new(void *tool_arguments)
{
  InterFeatArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file ...]", "Add intermediary "
                            "features between outside features in given "
                            "GFF3 file(s).");

  /* -outside */
  option = gt_option_new_string("outside", "set outside type",
                                arguments->outside_type, gt_ft_exon);
  gt_option_parser_add_option(op, option);

  /* -inter */
  option = gt_option_new_string("inter", "set intermediary type",
                                arguments->inter_type, gt_ft_intron);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  return op;
}

static int gt_interfeat_runner(int argc, const char **argv, int parsed_args,
                          void *tool_arguments, GtError *err)
{
  InterFeatArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream,
               *inter_feature_stream,
               *gff3_out_stream;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);

  /* create intermediary feature stream */
  inter_feature_stream =
    gt_inter_feature_stream_new(gff3_in_stream,
                                gt_str_get(arguments->outside_type),
                                gt_str_get(arguments->inter_type));

  /* create gff3 output stream */
  gff3_out_stream = gt_gff3_out_stream_new(inter_feature_stream,
                                           arguments->outfp);

  /* pull the features through the stream and free them afterwards */
  had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(inter_feature_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_interfeat(void)
{
  return gt_tool_new(gt_interfeat_arguments_new,
                     gt_interfeat_arguments_delete,
                     gt_interfeat_option_parser_new,
                     NULL,
                     gt_interfeat_runner);
}
