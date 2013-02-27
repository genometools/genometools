/*
  Copyright (c) 2008-2009, 2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/versionfunc.h"
#include "extended/bed_parser.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/bed_in_stream_api.h"
#include "tools/gt_bed_to_gff3.h"

typedef struct {
  GtStr *feature_type,
        *thick_feature_type,
        *block_type;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} BEDToGFF3Arguments;

static void *gt_bed_to_gff3_arguments_new(void)
{
  BEDToGFF3Arguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->feature_type = gt_str_new();
  arguments->thick_feature_type = gt_str_new();
  arguments->block_type = gt_str_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_bed_to_gff3_arguments_delete(void *tool_arguments)
{
  BEDToGFF3Arguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_delete(arguments->block_type);
  gt_str_delete(arguments->thick_feature_type);
  gt_str_delete(arguments->feature_type);
  gt_free(arguments);
}

static GtOptionParser* gt_bed_to_gff3_option_parser_new(void *tool_arguments)
{
  BEDToGFF3Arguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o;

  op = gt_option_parser_new("[BED_file]",
                            "Parse BED file and convert it to GFF3.");

  o = gt_option_new_string("featuretype", "Set type of parsed BED features",
                           arguments->feature_type, BED_FEATURE_TYPE);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_string("thicktype", "Set type of parsed thick BED features",
                           arguments->thick_feature_type,
                           BED_THICK_FEATURE_TYPE);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_string("blocktype", "Set type of parsed BED blocks",
                           arguments->block_type, BED_BLOCK_TYPE);
  gt_option_parser_add_option(op, o);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_max_args(op, 1);

  return op;
}

static int gt_bed_to_gff3_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args, void *tool_arguments,
                                 GtError *err)
{
  GtNodeStream *bed_in_stream = NULL, *gff3_out_stream = NULL;
  BEDToGFF3Arguments *arguments = tool_arguments;
  int had_err;

  gt_error_check(err);

  /* create a BED input stream */
  bed_in_stream = gt_bed_in_stream_new(argv[parsed_args]);
  gt_bed_in_stream_set_feature_type((GtBEDInStream*) bed_in_stream,
                                    gt_str_get(arguments->feature_type));
  gt_bed_in_stream_set_thick_feature_type((GtBEDInStream*) bed_in_stream,
                                          gt_str_get(arguments
                                                     ->thick_feature_type));
  gt_bed_in_stream_set_block_type((GtBEDInStream*) bed_in_stream,
                                  gt_str_get(arguments->block_type));

  /* create a GFF3 output stream */
  gff3_out_stream = gt_gff3_out_stream_new(bed_in_stream, arguments->outfp);

  /* pull the features through the stream and free them afterwards */
  had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(bed_in_stream);

  return had_err;
}

GtTool* gt_bed_to_gff3(void)
{
  return gt_tool_new(gt_bed_to_gff3_arguments_new,
                     gt_bed_to_gff3_arguments_delete,
                     gt_bed_to_gff3_option_parser_new,
                     NULL,
                     gt_bed_to_gff3_runner);
}
