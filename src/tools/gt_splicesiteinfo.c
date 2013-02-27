/*
  Copyright (c) 2007-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/versionfunc.h"
#include "core/warning_api.h"
#include "extended/add_introns_stream_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/splice_site_info_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "tools/gt_splicesiteinfo.h"

typedef struct {
  bool addintrons;
  GtSeqid2FileInfo *s2fi;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} SpliceSiteInfoArguments;

static void* gt_splicesiteinfo_arguments_new(void)
{
  SpliceSiteInfoArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_splicesiteinfo_arguments_delete(void *tool_arguments)
{
  SpliceSiteInfoArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_splicesiteinfo_option_parser_new(void *tool_arguments)
{
  SpliceSiteInfoArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [GFF3_file ...]", "Show information "
                            "about splice sites given in GFF3 files.");

  /* -seqfile, -matchdesc, -usedesc and -regionmapping */
  gt_seqid2file_register_options(op, arguments->s2fi);

  /* -addintrons */
  option = gt_option_new_bool("addintrons", "add intron features between "
                              "existing exon features\n(before computing the "
                              "information to be shown)",
                              &arguments->addintrons, false);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);

  return op;
}

static int gt_splicesiteinfo_runner(int argc, const char **argv,
                                    int parsed_args, void *tool_arguments,
                                    GtError *err)
{
  SpliceSiteInfoArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream = NULL,
               *add_introns_stream = NULL,
               *splice_site_info_stream = NULL;
  GtRegionMapping *region_mapping;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (!had_err) {
    /* create gff3 input stream */
    gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                    argv + parsed_args);

    /* create region mapping */
    region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!region_mapping)
      had_err = -1;
  }

  if (!had_err) {
    /* create addintrons stream (if necessary) */
    if (arguments->addintrons)
      add_introns_stream = gt_add_introns_stream_new(gff3_in_stream);

    /* create extract feature stream */
    splice_site_info_stream = gt_splice_site_info_stream_new(
                                                          arguments->addintrons
                                                          ? add_introns_stream
                                                          : gff3_in_stream,
                                                          region_mapping);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(splice_site_info_stream, err);
  }

  if (!had_err) {
    if (!gt_splice_site_info_stream_show(splice_site_info_stream,
                                         arguments->outfp)) {
      gt_warning("input file(s) contained no intron, use option -addintrons to "
                 "add introns automatically");
    }
  }

  /* free */
  gt_node_stream_delete(splice_site_info_stream);
  gt_node_stream_delete(add_introns_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_splicesiteinfo(void)
{
  return gt_tool_new(gt_splicesiteinfo_arguments_new,
                     gt_splicesiteinfo_arguments_delete,
                     gt_splicesiteinfo_option_parser_new,
                     NULL,
                     gt_splicesiteinfo_runner);
}
