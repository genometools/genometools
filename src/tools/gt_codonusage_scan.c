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
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "extended/codon_usage_visitor.h"
#include "extended/codon_usage_scan_stream.h"
#include "extended/visitor_stream.h"
#include "tools/gt_codonusage_scan.h"

typedef struct {
  GtSeqid2FileInfo *s2fi;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  GtStr *freqfile;
  unsigned int windowsize;
  double coding_threshold;
  bool all;
} GtCodonusageScanArguments;

static void* gt_codonusage_scan_arguments_new(void)
{
  GtCodonusageScanArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  arguments->freqfile = gt_str_new();
  return arguments;
}

static void gt_codonusage_scan_arguments_delete(void *tool_arguments)
{
  GtCodonusageScanArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_str_delete(arguments->freqfile);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_codonusage_scan_option_parser_new(void
                                                                *tool_arguments)
{
  GtCodonusageScanArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [GFF3_file]",
                            "Assess coding potential of CDS features "
                            "by codon usage.");

  option = gt_option_new_filename("freqfile", "file with codon frequencies",
                                  arguments->freqfile);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  option = gt_option_new_uint("windowsize", "window size for codon "
                                            "usage analysis",
                              &arguments->windowsize, 100);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_double("threshold", "threshold for coding "
                                             "classification",
                                &arguments->coding_threshold, 4.0);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("all", "show all ORFs (i.e. do not filter)",
                                &arguments->all, false);
  gt_option_parser_add_option(op, option);

  /* -seqfile, -matchdesc, -usedesc and -regionmapping */
  gt_seqid2file_register_options(op, arguments->s2fi);

  /* output file options */
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  gt_option_parser_set_max_args(op, 1);

  return op;
}

static int gt_codonusage_scan_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args, void *tool_arguments,
                                 GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL,
               *codonusage_scan_stream1 = NULL,
               *codonusage_scan_stream2 = NULL,
               *gff3_out_stream = NULL,
               *last_stream = NULL;
  GtCodonusageScanArguments *arguments = tool_arguments;
  GtRegionMapping *region_mapping = NULL;
  GtNodeVisitor *nv = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (!had_err) {
    /* create gff3 input stream */
    last_stream = gff3_in_stream =
                                gt_gff3_in_stream_new_sorted(argv[parsed_args]);

    /* create region mapping */
    region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!region_mapping)
      had_err = -1;
  }

  if (!had_err) {
    nv = gt_codon_usage_visitor_new(region_mapping,
                                    gt_str_get(arguments->freqfile),
                                    arguments->windowsize,
                                    arguments->coding_threshold,
                                    err);
    if (!nv)
      had_err = -1;
  }
  if (!had_err) {
    last_stream = codonusage_scan_stream1 = gt_visitor_stream_new(last_stream,
                                                                  nv);
    if (!arguments->all)
      last_stream = codonusage_scan_stream2 =
                                    gt_codon_usage_scan_stream_new(last_stream);

    last_stream = gff3_out_stream = gt_gff3_out_stream_new(last_stream,
                                                           arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(last_stream, err);
  }

  /* free */
  gt_node_stream_delete(codonusage_scan_stream1);
  gt_node_stream_delete(codonusage_scan_stream2);
  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(gff3_out_stream);
  gt_region_mapping_delete(region_mapping);

  return had_err;
}

GtTool* gt_codonusage_scan(void)
{
  return gt_tool_new(gt_codonusage_scan_arguments_new,
                     gt_codonusage_scan_arguments_delete,
                     gt_codonusage_scan_option_parser_new,
                     NULL,
                     gt_codonusage_scan_runner);
}
