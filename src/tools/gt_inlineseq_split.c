/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd

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

#include "core/file_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/node_stream_api.h"
#include "extended/gff3_in_stream_api.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/seqid2file.h"
#include "extended/sequence_node_out_stream.h"
#include "extended/sequence_node_out_visitor.h"

#include "tools/gt_inlineseq_split.h"

typedef struct {
  GtStr *seqoutfile,
        *gffoutfile;
} GtInlineseqSplitArguments;

static void* gt_inlineseq_split_arguments_new(void)
{
  GtInlineseqSplitArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->seqoutfile = gt_str_new();
  arguments->gffoutfile = gt_str_new();
  return arguments;
}

static void gt_inlineseq_split_arguments_delete(void *tool_arguments)
{
  GtInlineseqSplitArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->seqoutfile);
    gt_str_delete(arguments->gffoutfile);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_inlineseq_split_option_parser_new(void *tool_arguments)
{
  GtInlineseqSplitArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption GT_UNUSED *option, *seqopt, *gff3opt;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("(-seqfile <seqfile> | -gff3file <gff3file>) "
                            "[GFF3_file]",
                            "Split GFF3 annotations with inline sequences into "
                            "separate files.");


  /* -seqfile */
  seqopt = gt_option_new_string("seqfile", "output file for sequences as FASTA",
                                arguments->seqoutfile, NULL);
  gt_option_parser_add_option(op, seqopt);

  /* -gff3file */
  gff3opt = gt_option_new_string("gff3file", "output file for annotations "
                                             "as GFF3",
                                 arguments->gffoutfile, NULL);
  gt_option_parser_add_option(op, gff3opt);

  gt_option_is_mandatory_either(seqopt, gff3opt);

  return op;
}

static int gt_inlineseq_split_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GtError *err)
{
  GtInlineseqSplitArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream = NULL,
               *gff3_out_stream = NULL,
               *split_stream = NULL,
               *last_stream = NULL;
  GtFile *seq_out_file = NULL,
         *gff3_out_file = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);


  if (gt_str_length(arguments->seqoutfile) > 0) {
    seq_out_file = gt_file_new(gt_str_get(arguments->seqoutfile), "w+", err);
    if (!seq_out_file)
      had_err = -1;
  }

  if (!had_err && gt_str_length(arguments->gffoutfile) > 0) {
    gff3_out_file = gt_file_new(gt_str_get(arguments->gffoutfile), "w+", err);
    if (!gff3_out_file)
      had_err = -1;
  }

  if (!had_err) {
    last_stream = gff3_in_stream = gt_gff3_in_stream_new_unsorted(
                                                            argc - parsed_args,
                                                            argv + parsed_args);
    gt_assert(gff3_in_stream);
  }

  if (!had_err) {
    last_stream = split_stream = gt_sequence_node_out_stream_new(last_stream,
                                                                 seq_out_file,
                                                                 err);
    gt_assert(split_stream);
  }

  if (!had_err) {
    last_stream = gff3_out_stream = gt_gff3_out_stream_new(last_stream,
                                                           gff3_out_file);
    had_err = gt_node_stream_pull(last_stream, err);
  }

  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(split_stream);
  gt_file_delete(seq_out_file);
  gt_file_delete(gff3_out_file);

  return had_err;
}

GtTool* gt_inlineseq_split(void)
{
  return gt_tool_new(gt_inlineseq_split_arguments_new,
                     gt_inlineseq_split_arguments_delete,
                     gt_inlineseq_split_option_parser_new,
                     NULL,
                     gt_inlineseq_split_runner);
}
