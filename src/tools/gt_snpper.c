/*
  Copyright (c) 2012-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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
#include "core/trans_table_api.h"
#include "core/unused_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "extended/snp_annotator_stream.h"
#include "tools/gt_snpper.h"

typedef struct {
  GtSeqid2FileInfo *s2fi;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  unsigned int ttable;
  GtStr *desc;
} GtSnpperArguments;

static void* gt_snpper_arguments_new(void)
{
  GtStrArray *descs;
  unsigned long i;
  GtSnpperArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->ofi = gt_output_file_info_new();
  arguments->desc = gt_str_new_cstr("NCBI translation table number, "
                                    "choose from:\n");
  descs = gt_trans_table_get_scheme_descriptions();
  for (i = 0; i < gt_str_array_size(descs); i++) {
    gt_str_append_cstr(arguments->desc, gt_str_array_get(descs, i));
    if (i != gt_str_array_size(descs)-1)
      gt_str_append_cstr(arguments->desc, "\n");
  }
  gt_str_array_delete(descs);
  return arguments;
}

static void gt_snpper_arguments_delete(void *tool_arguments)
{
  GtSnpperArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_str_delete(arguments->desc);
  gt_output_file_info_delete(arguments->ofi);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_free(arguments);
}

static GtOptionParser* gt_snpper_option_parser_new(void *tool_arguments)
{
  GtSnpperArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] GFF3_file [GVF_file]",
                            "Annotates SNPs according to their effect on "
                            "the genome as given by a genomic annotation.");

  /* -trans_table */
  option = gt_option_new_uint("trans_table", gt_str_get(arguments->desc),
                              &arguments->ttable, 1);
  gt_option_parser_add_option(op, option);

  /* -seqfile, -encseq, -matchdesc, -usedesc and -regionmapping */
  gt_seqid2file_register_options(op, arguments->s2fi);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_min_max_args(op, 1, 2);

  return op;
}

static int gt_snpper_runner(GT_UNUSED int argc, const char **argv,
                                 int parsed_args, void *tool_arguments,
                                 GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL,
               *gvf_in_stream = NULL,
               *gvf_out_stream = NULL,
               *snp_annotator_stream = NULL;
  GtSnpperArguments *arguments = tool_arguments;
  GtRegionMapping *region_mapping = NULL;
  GtTransTable *tt = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (!(tt = gt_trans_table_new(arguments->ttable, err)))
    had_err = -1;

  if (!had_err) {
    /* create GFF3 input stream */
    gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);
    if (!gff3_in_stream)
      had_err = -1;
  }

  if (!had_err) {
    /* create GVF input stream */
    gvf_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args+1]);
    if (!gvf_in_stream)
      had_err = -1;
  }

  if (!had_err) {
    /* create region mapping */
    region_mapping = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!region_mapping)
      had_err = -1;
  }

  if (!had_err) {
    gt_assert(tt);
    snp_annotator_stream = gt_snp_annotator_stream_new(gvf_in_stream,
                                                       gff3_in_stream,
                                                       tt,
                                                       region_mapping);
    if (!snp_annotator_stream)
      had_err = -1;
  }

  if (!had_err) {
    gvf_out_stream = gt_gff3_out_stream_new(snp_annotator_stream,
                                            arguments->outfp);
    if (!gvf_out_stream)
      had_err = -1;
  }

  if (!had_err) {
    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(gvf_out_stream, err);
  }

  /* free */
  gt_node_stream_delete(snp_annotator_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(gvf_in_stream);
  gt_node_stream_delete(gvf_out_stream);
  gt_trans_table_delete(tt);
  gt_region_mapping_delete(region_mapping);

  return had_err;
}

GtTool* gt_snpper(void)
{
  return gt_tool_new(gt_snpper_arguments_new,
                     gt_snpper_arguments_delete,
                     gt_snpper_option_parser_new,
                     NULL,
                     gt_snpper_runner);
}
