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

#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gtf_out_stream_api.h"
#include "tools/gt_gff3_to_gtf.h"

typedef struct {
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GFF3ToGTFArguments;

static void* gt_gff3_to_gtf_arguments_new(void)
{
  GFF3ToGTFArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_gff3_to_gtf_arguments_delete(void *tool_arguments)
{
  GFF3ToGTFArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_gff3_to_gtf_option_parser_new(void *tool_arguments)
{
  GFF3ToGTFArguments *arguments = tool_arguments;
  GtOptionParser *op;
  op = gt_option_parser_new("[GFF3_file ...]",
                            "Parse GFF3 file(s) and show them as GTF2.2.");
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);
  return op;
}

static int gt_gff3_to_gtf_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GFF3ToGTFArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream = NULL, *gtf_out_stream = NULL;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);

  /* create a gtf output stream */
  gtf_out_stream = gt_gtf_out_stream_new(gff3_in_stream, arguments->outfp);

  /* pull the features through the stream and free them afterwards */
  had_err = gt_node_stream_pull(gtf_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(gtf_out_stream);

  return had_err;
}

GtTool* gt_gff3_to_gtf(void)
{
  return gt_tool_new(gt_gff3_to_gtf_arguments_new,
                     gt_gff3_to_gtf_arguments_delete,
                     gt_gff3_to_gtf_option_parser_new,
                     NULL,
                     gt_gff3_to_gtf_runner);
}
