/*
  Copyright (c) 2010-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "extended/md5_to_id_stream_api.h"
#include "tools/gt_md5_to_id.h"

typedef struct {
  bool verbose;
  GtSeqid2FileInfo *s2fi;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} MD5ToSeqidsArguments;

static void *gt_md5_to_id_arguments_new(void)
{
  MD5ToSeqidsArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_md5_to_id_arguments_delete(void *tool_arguments)
{
  MD5ToSeqidsArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_md5_to_id_option_parser_new(void *tool_arguments)
{
  MD5ToSeqidsArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                            "Change MD5 fingerprints used as sequence IDs in "
                            "given GFF3 files to ``regular'' ones.");

  /* -seqfile, -matchdesc, -usedesc and -regionmapping */
  gt_seqid2file_register_options_ext(op, arguments->s2fi, false, true);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  return op;
}

static int gt_md5_to_id_runner(GT_UNUSED int argc, const char **argv,
                                   int parsed_args, void *tool_arguments,
                                   GtError *err)
{
  GtNodeStream *gff3_in_stream, *md5_to_id_stream = NULL,
               *gff3_out_stream = NULL;
  MD5ToSeqidsArguments *arguments = tool_arguments;
  GtRegionMapping *region_mapping = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  if (arguments->verbose && arguments->outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create region mapping */
  if (gt_seqid2file_option_used(arguments->s2fi)) {
    region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!region_mapping)
      had_err = -1;
  }

  if (!had_err) {
    /* create seqid to md5 stream */
    md5_to_id_stream = gt_md5_to_id_stream_new(gff3_in_stream, region_mapping);

    /* create gff3 output stream */
    gff3_out_stream = gt_gff3_out_stream_new(md5_to_id_stream,
                                             arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(gff3_out_stream, err);
  }

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(md5_to_id_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool *gt_md5_to_id(void)
{
  return gt_tool_new(gt_md5_to_id_arguments_new,
                     gt_md5_to_id_arguments_delete,
                     gt_md5_to_id_option_parser_new,
                     NULL,
                     gt_md5_to_id_runner);
}
