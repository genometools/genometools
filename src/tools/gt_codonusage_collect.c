/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>

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
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "extended/codon_usage_collector.h"
#include "extended/visitor_stream.h"
#include "tools/gt_codonusage_collect.h"

typedef struct {
  GtSeqid2FileInfo *s2fi;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GtCodonusageCollectArguments;

static void* gt_codonusage_collect_arguments_new(void)
{
  GtCodonusageCollectArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_codonusage_collect_arguments_delete(void *tool_arguments)
{
  GtCodonusageCollectArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_codonusage_collect_option_parser_new(void
                                                                *tool_arguments)
{
  GtCodonusageCollectArguments *arguments = tool_arguments;
  GtOptionParser *op;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [GFF3_file]",
                            "Calculates codon frequencies for "
                            "coding sequences.");

  /* -seqfile, -matchdesc, -usedesc and -regionmapping */
  gt_seqid2file_register_options(op, arguments->s2fi);

  /* output file options */
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  gt_option_parser_set_max_args(op, 1);

  return op;
}

static int gt_codonusage_collect_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args, void *tool_arguments,
                                 GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL, *codonusage_collect_stream = NULL;
  GtCodonusageCollectArguments *arguments = tool_arguments;
  GtRegionMapping *region_mapping;
  GtNodeVisitor *nv;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (!had_err) {
    /* create gff3 input stream */
    gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);

    /* create region mapping */
    region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!region_mapping)
      had_err = -1;
  }

  if (!had_err) {
    nv = gt_codon_usage_collector_new(region_mapping);
    gt_assert(nv);
    codonusage_collect_stream = gt_visitor_stream_new(gff3_in_stream, nv);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(codonusage_collect_stream, err);
  }

  if (!had_err)
    gt_codon_usage_collector_output_usage_table((GtCodonUsageCollector*) nv,
                                                arguments->outfp);

  /* free */
  gt_node_stream_delete(codonusage_collect_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_region_mapping_delete(region_mapping);

  return had_err;
}

GtTool* gt_codonusage_collect(void)
{
  return gt_tool_new(gt_codonusage_collect_arguments_new,
                     gt_codonusage_collect_arguments_delete,
                     gt_codonusage_collect_option_parser_new,
                     NULL,
                     gt_codonusage_collect_runner);
}
