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
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/seqid2file.h"
#include "extended/sequence_node_add_stream.h"
#include "extended/visitor_stream.h"
#include "tools/gt_inlineseq_add.h"

typedef struct {
  GtSeqid2FileInfo *s2fi;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} InlineseqAddArguments;

static void *gt_inlineseq_add_arguments_new(void)
{
  InlineseqAddArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  arguments->outfp = NULL;
  return arguments;
}

static void gt_inlineseq_add_arguments_delete(void *tool_arguments)
{
  InlineseqAddArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_inlineseq_add_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GT_UNUSED GtOption *option;
  InlineseqAddArguments *arguments = tool_arguments;

  /* init */
  op = gt_option_parser_new("[options] [GFF3_file ...]",
                            "Adds inline sequences from external source to "
                            "GFF3 input.");

  gt_seqid2file_register_options(op, arguments->s2fi);
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  return op;
}

static int gt_inlineseq_add_runner(int argc, const char **argv, int parsed_args,
                               void *tool_arguments, GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL,
               *add_stream = NULL,
               *gff3_out_stream = NULL,
               *last_stream = NULL;
  GtRegionMapping *rm = NULL;
  InlineseqAddArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);

  /* add region mapping if given */
  if (gt_seqid2file_option_used(arguments->s2fi)) {
    rm = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!rm)
      had_err = -1;
  }

  last_stream = gff3_in_stream = gt_gff3_in_stream_new_unsorted(
                                                            argc - parsed_args,
                                                            argv + parsed_args);
  gt_assert(gff3_in_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream*) gff3_in_stream);

  last_stream = add_stream = gt_sequence_node_add_stream_new(last_stream, rm,
                                                             err);
  if (!add_stream) {
    had_err = -1;
  }

  if (!had_err) {
    last_stream = gff3_out_stream = gt_gff3_out_stream_new(last_stream,
                                                           arguments->outfp);
  }

  if (!had_err)
    had_err = gt_node_stream_pull(last_stream, err);

  /* free */
  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(add_stream);
  gt_node_stream_delete(gff3_out_stream);
  gt_region_mapping_delete(rm);

  return had_err;
}

GtTool* gt_inlineseq_add(void)
{
  return gt_tool_new(gt_inlineseq_add_arguments_new,
                     gt_inlineseq_add_arguments_delete,
                     gt_inlineseq_add_option_parser_new,
                     NULL,
                     gt_inlineseq_add_runner);
}
