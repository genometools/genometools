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
#include "extended/cds_stream_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "tools/gt_cds.h"

#define GT_CDS_SOURCE_TAG "gt cds"

typedef struct {
  unsigned int minorflen;
  bool start_codon,
       final_stop_codon,
       generic_start_codons,
       verbose;
  GtSeqid2FileInfo *s2fi;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} CDSArguments;

static void* gt_cds_arguments_new(void)
{
  CDSArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_cds_arguments_delete(void *tool_arguments)
{
  CDSArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_cds_option_parser_new(void *tool_arguments)
{
  CDSArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] [GFF3_file]",
                            "Add CDS (coding sequence) features to exon "
                            "features given in GFF3 file.");

  /* -minorflen */
  option = gt_option_new_uint_min("minorflen", "set the minimum length an open "
                                  "reading frame (ORF) must have to be added "
                                  "as a CDS feature (measured in amino acids)",
                                  &arguments->minorflen, 64, 1);
  gt_option_parser_add_option(op, option);

  /* -startcodon */
  option = gt_option_new_bool("startcodon", "require than an ORF must begin "
                              "with a start codon", &arguments->start_codon,
                              false);
  gt_option_parser_add_option(op, option);

  /* -finalstopcodon */
  option = gt_option_new_bool("finalstopcodon", "require that the final ORF "
                              "must end with a stop codon",
                              &arguments->final_stop_codon, false);
  gt_option_parser_add_option(op, option);

  /* -genericstartcodons */
  option = gt_option_new_bool("genericstartcodons", "use generic start codons",
                              &arguments->generic_start_codons, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -seqfile, -matchdesc, -usedesc and -regionmapping */
  gt_seqid2file_register_options(op, arguments->s2fi);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_max_args(op, 1);

  return op;
}

static int gt_cds_runner(GT_UNUSED int argc, const char **argv, int parsed_args,
                         void *tool_arguments, GtError *err)
{
  GtNodeStream *gff3_in_stream, *cds_stream = NULL, *gff3_out_stream = NULL;
  CDSArguments *arguments = tool_arguments;
  GtRegionMapping *region_mapping;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* create gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);
  if (arguments->verbose && arguments->outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create region mapping */
  region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
  if (!region_mapping)
    had_err = -1;

  if (!had_err) {
    /* create CDS stream */
    cds_stream = gt_cds_stream_new(gff3_in_stream, region_mapping,
                                   arguments->minorflen, GT_CDS_SOURCE_TAG,
                                   arguments->start_codon,
                                   arguments->final_stop_codon,
                                   arguments->generic_start_codons);

    /* create gff3 output stream */
    gff3_out_stream = gt_gff3_out_stream_new(cds_stream, arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(gff3_out_stream, err);
  }

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(cds_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_cds(void)
{
  return gt_tool_new(gt_cds_arguments_new,
                     gt_cds_arguments_delete,
                     gt_cds_option_parser_new,
                     NULL,
                     gt_cds_runner);
}
