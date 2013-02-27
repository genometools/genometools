/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "extended/extract_feature_stream_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "tools/gt_extractfeat.h"

typedef struct {
  bool join,
       translate,
       seqid,
       target,
       verbose;
  GtStr *type;
  GtSeqid2FileInfo *s2fi;
  unsigned long width;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GtExtractFeatArguments;

static void* gt_extractfeat_arguments_new(void)
{
  GtExtractFeatArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->type = gt_str_new();
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_extractfeat_arguments_delete(void *tool_arguments)
{
  GtExtractFeatArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_str_delete(arguments->type);
  gt_free(arguments);
}

static GtOptionParser* gt_extractfeat_option_parser_new(void *tool_arguments)
{
  GtExtractFeatArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [GFF3_file]",
                            "Extract features given in GFF3 file from "
                            "sequence file.");

  /* -type */
  option = gt_option_new_string("type", "set type of features to extract",
                                arguments->type, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* -join */
  option = gt_option_new_bool("join", "join feature sequences in the same "
                              "subgraph into a single one", &arguments->join,
                              false);
  gt_option_parser_add_option(op, option);

  /* -translate */
  option = gt_option_new_bool("translate", "translate the features (of a DNA "
                              "sequence) into protein", &arguments->translate,
                              false);
  gt_option_parser_add_option(op, option);

  /* -seqid */
  option = gt_option_new_bool("seqid", "add sequence ID of extracted features "
                              "to FASTA descriptions", &arguments->seqid,
                              false);
  gt_option_parser_add_option(op, option);

  /* -target */
  option = gt_option_new_bool("target", "add target ID(s) of extracted "
                              "features to FASTA descriptions",
                              &arguments->target, false);
  gt_option_parser_add_option(op, option);

  /* -seqfile, -matchdesc, -usedesc and -regionmapping */
  gt_seqid2file_register_options(op, arguments->s2fi);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* -width */
  option = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_max_args(op, 1);

  return op;
}

static int gt_extractfeat_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args, void *tool_arguments,
                                 GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL, *extract_feature_stream = NULL;
  GtExtractFeatArguments *arguments = tool_arguments;
  GtRegionMapping *region_mapping;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (!had_err) {
    /* create gff3 input stream */
    gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);
    if (arguments->verbose)
      gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

    /* create region mapping */
    region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!region_mapping)
      had_err = -1;
  }

  if (!had_err) {
    /* create extract feature stream */
    extract_feature_stream =
      gt_extract_feature_stream_new(gff3_in_stream,
                                    region_mapping,
                                    gt_str_get(arguments->type),
                                    arguments->join,
                                    arguments->translate,
                                    arguments->seqid,
                                    arguments->target,
                                    arguments->width,
                                    arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(extract_feature_stream, err);
  }

  /* free */
  gt_node_stream_delete(extract_feature_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_extractfeat(void)
{
  return gt_tool_new(gt_extractfeat_arguments_new,
                     gt_extractfeat_arguments_delete,
                     gt_extractfeat_option_parser_new,
                     NULL,
                     gt_extractfeat_runner);
}
