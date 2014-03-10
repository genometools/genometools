/*
  Copyright (c) 2009 Gordon Gremme <gordon@gremme.org>

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
#include "extended/dup_feature_stream_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "tools/gt_interfeat.h"

typedef struct {
  GtStr *dest_type,
        *source_type;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} InterFeatArguments;

static void* gt_interfeat_arguments_new(void)
{
  InterFeatArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->dest_type = gt_str_new();
  arguments->source_type = gt_str_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_interfeat_arguments_delete(void *tool_arguments)
{
  InterFeatArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_delete(arguments->source_type);
  gt_str_delete(arguments->dest_type);
  gt_free(arguments);
}

static GtOptionParser* gt_interfeat_option_parser_new(void *tool_arguments)
{
  InterFeatArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file ...]", "Duplicate "
                            "internal feature nodes in given GFF3 files.");

  /* -outside */
  option = gt_option_new_string("dest", "set destination type",
                                arguments->dest_type, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* -inter */
  option = gt_option_new_string("source", "set source type",
                                arguments->source_type, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  return op;
}

static int gt_interfeat_runner(int argc, const char **argv, int parsed_args,
                          void *tool_arguments, GtError *err)
{
  InterFeatArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream,
               *dup_feature_stream,
               *gff3_out_stream;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);

  /* create intermediary feature stream */
  dup_feature_stream =
    gt_dup_feature_stream_new(gff3_in_stream, gt_str_get(arguments->dest_type),
                              gt_str_get(arguments->source_type));

  /* create gff3 output stream */
  gff3_out_stream = gt_gff3_out_stream_new(dup_feature_stream,
                                           arguments->outfp);

  /* pull the features through the stream and free them afterwards */
  had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(dup_feature_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_dupfeat(void)
{
  return gt_tool_new(gt_interfeat_arguments_new,
                     gt_interfeat_arguments_delete,
                     gt_interfeat_option_parser_new,
                     NULL,
                     gt_interfeat_runner);
}
