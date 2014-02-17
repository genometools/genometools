/*
  Copyright (c) 2006-2011, 2013 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008       Center for Bioinformatics, University of Hamburg

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
#include "core/versionfunc.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/merge_stream_api.h"
#include "tools/gt_merge.h"

typedef struct {
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  bool retainids;
} MergeArguments;

static void* gt_merge_arguments_new(void)
{
  MergeArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_merge_arguments_delete(void *tool_arguments)
{
  MergeArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_merge_option_parser_new(void *tool_arguments)
{
  MergeArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);
  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                         "Merge sorted GFF3 files in sorted fashion.");

  /* -retainids */
  option = gt_option_new_bool("retainids",
                              "when available, use the original IDs provided "
                              "in the source file\n"
                              "(memory consumption is proportional to the "
                              "input file size(s))", &arguments->retainids,
                              false);
  gt_option_parser_add_option(op, option);

  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);
  return op;
}

static int gt_merge_runner(int argc, const char **argv, int parsed_args,
                           void *tool_arguments, GtError *err)
{
  MergeArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream,
                *merge_stream,
                *gff3_out_stream;
  GtArray *genome_streams;
  GtUword i;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  /* alloc */
  genome_streams = gt_array_new(sizeof (GtNodeStream*));

  /* XXX: check for multiple specification of '-' */

  /* create an gff3 input stream for each given file */
  if (parsed_args < argc) {
    /* we got files to open */
    for (i = parsed_args; i < argc; i++) {
      gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[i]);
      gt_array_add(genome_streams, gff3_in_stream);
    }
   }
   else {
     /* use stdin */
     gff3_in_stream = gt_gff3_in_stream_new_sorted(NULL);
     gt_array_add(genome_streams, gff3_in_stream);
   }

  /* create a merge stream */
  merge_stream = gt_merge_stream_new(genome_streams);
  gt_assert(merge_stream);

  /* create a gff3 output stream */
  gff3_out_stream = gt_gff3_out_stream_new(merge_stream, arguments->outfp);
  if (arguments->retainids)
    gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream*) gff3_out_stream);

  /* pull the features through the stream and free them afterwards */
  had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(merge_stream);
  for (i = 0; i < gt_array_size(genome_streams); i++)
    gt_node_stream_delete(*(GtNodeStream**) gt_array_get(genome_streams, i));
  gt_array_delete(genome_streams);

  return had_err;
}

GtTool* gt_merge(void)
{
  return gt_tool_new(gt_merge_arguments_new,
                     gt_merge_arguments_delete,
                     gt_merge_option_parser_new,
                     NULL,
                     gt_merge_runner);
}
