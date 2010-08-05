/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/option.h"
#include "core/outputfile.h"
#include "core/unused_api.h"
#include "extended/cds_stream.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/seqid2file.h"
#include "extended/seqids_to_md5_stream.h"
#include "tools/gt_seqids_to_md5.h"

typedef struct {
  bool verbose,
       usedesc;
  GtStr *seqfile,
        *region_mapping;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} SeqidsToMD5Arguments;

static void *gt_seqids_to_md5_arguments_new(void)
{
  SeqidsToMD5Arguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->seqfile = gt_str_new();
  arguments->region_mapping = gt_str_new();
  arguments->ofi = gt_outputfileinfo_new();
  return arguments;
}

static void gt_seqids_to_md5_arguments_delete(void *tool_arguments)
{
  SeqidsToMD5Arguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_str_delete(arguments->region_mapping);
  gt_str_delete(arguments->seqfile);
  gt_free(arguments);
}

static GtOptionParser* gt_seqids_to_md5_option_parser_new(void *tool_arguments)
{
  SeqidsToMD5Arguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                            "Change sequence IDs in given GFF3 files to MD5 "
                            "fingerprints of the corresponding sequences.");

  /* -seqfile, -usedesc and -regionmapping */
  gt_seqid2file_options(op, arguments->seqfile, &arguments->usedesc,
                        arguments->region_mapping);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  return op;
}

static int gt_seqids_to_md5_runner(GT_UNUSED int argc, const char **argv,
                                   int parsed_args, void *tool_arguments,
                                   GtError *err)
{
  GtNodeStream *gff3_in_stream, *seqids_to_md5_stream = NULL,
               *gff3_out_stream = NULL;
  SeqidsToMD5Arguments *arguments = tool_arguments;
  GtRegionMapping *region_mapping;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  if (arguments->verbose && arguments->outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create region mapping */
  region_mapping = gt_seqid2file_regionmapping_new(arguments->seqfile,
                                                   arguments->usedesc,
                                                   arguments->region_mapping,
                                                   err);
  if (!region_mapping)
    had_err = -1;

  if (!had_err) {
    /* create seqid to md5 stream */
    seqids_to_md5_stream = gt_seqids_to_md5_stream_new(gff3_in_stream,
                                                       region_mapping);

    /* create gff3 output stream */
    gff3_out_stream = gt_gff3_out_stream_new(seqids_to_md5_stream,
                                             arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(gff3_out_stream, err);
  }

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(seqids_to_md5_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool *gt_seqids_to_md5(void)
{
  return gt_tool_new(gt_seqids_to_md5_arguments_new,
                     gt_seqids_to_md5_arguments_delete,
                     gt_seqids_to_md5_option_parser_new,
                     NULL,
                     gt_seqids_to_md5_runner);
}
