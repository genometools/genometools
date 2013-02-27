/*
  Copyright (c) 2005-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gff3_output.h"
#include "extended/gff3_parser.h"
#include "extended/gff3_visitor.h"
#include "extended/gtdatahelp.h"
#include "extended/select_stream.h"
#include "extended/targetbest_select_stream.h"
#include "tools/gt_select.h"

#define GT_STRAND_OPT  "strand"
#define TARGETGT_STRAND_OPT  "targetstrand"

typedef struct {
  bool verbose,
       has_CDS,
       targetbest;
  GtStr *seqid,
        *source,
        *gt_strand_char,
        *targetgt_strand_char;
  GtRange contain_range,
          overlap_range;
  GtStrand strand,
           targetstrand;
  unsigned long max_gene_length,
                max_gene_num,
                feature_num;
  double min_gene_score,
         max_gene_score,
         min_average_splice_site_prob,
         single_intron_factor;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  GtStrArray  *filter_files;
  GtStr *filter_logic;
  GtStr *dropped_file;
} SelectArguments;

static void* gt_select_arguments_new(void)
{
  SelectArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->seqid = gt_str_new();
  arguments->source = gt_str_new();
  arguments->gt_strand_char = gt_str_new();
  arguments->strand = GT_NUM_OF_STRAND_TYPES;
  arguments->targetgt_strand_char = gt_str_new();
  arguments->targetstrand = GT_NUM_OF_STRAND_TYPES;
  arguments->ofi = gt_output_file_info_new();
  arguments->filter_files = gt_str_array_new();
  arguments->filter_logic = gt_str_new();
  arguments->dropped_file = gt_str_new();
  return arguments;
}

static void gt_select_arguments_delete(void *tool_arguments)
{
  SelectArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_delete(arguments->targetgt_strand_char);
  gt_str_delete(arguments->gt_strand_char);
  gt_str_delete(arguments->source);
  gt_str_delete(arguments->seqid);
  gt_str_array_delete(arguments->filter_files);
  gt_str_delete(arguments->filter_logic);
  gt_str_delete(arguments->dropped_file);
  gt_free(arguments);
}

static GtOptionParser* gt_select_option_parser_new(void *tool_arguments)
{
  SelectArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *contain_option, *overlap_option, *minaveragessp_option,
           *singleintron_option, *optiondroppedfile;
  gt_assert(arguments);

  static const char *filter_logic[] = {
    "AND",
    "OR",
    NULL
  };

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                            "Select certain features (specified by the used "
                            "options) from given GFF3 file(s).");

  /* -seqid */
  option = gt_option_new_string("seqid", "select feature with the given "
                                "sequence ID (all comments are selected). ",
                                arguments->seqid, NULL);
  gt_option_parser_add_option(op, option);

  /* -source */
  option = gt_option_new_string("source", "select feature with the given "
                                "source (the source is column 2 in regular "
                                "GFF3 lines)" , arguments->source, NULL);
  gt_option_parser_add_option(op, option);

  /* -contain */
  contain_option = gt_option_new_range("contain", "select all features which "
                                       "are contained in the given range",
                                       &arguments->contain_range, NULL);
  gt_option_parser_add_option(op, contain_option);

  /* -overlap */
  overlap_option = gt_option_new_range("overlap", "select all features which "
                                       "do overlap with the given range",
                                       &arguments->overlap_range, NULL);
  gt_option_parser_add_option(op, overlap_option);

  /* -strand */
  option = gt_option_new_string(GT_STRAND_OPT, "select all top-level features"
                                "(i.e., features without parents) whose strand "
                                "equals the given one (must be one of '"
                                GT_STRAND_CHARS"')", arguments->gt_strand_char,
                                NULL);
  gt_option_parser_add_option(op, option);

  /* -targetstrand */
  option = gt_option_new_string(TARGETGT_STRAND_OPT, "select all top-level "
                                "features (i.e., features without parents) "
                                "which have exactly one target attribute whose "
                                "strand equals the given one (must be one of '"
                                GT_STRAND_CHARS"')",
                                arguments->targetgt_strand_char, NULL);
  gt_option_parser_add_option(op, option);

  /* -targetbest */
  option = gt_option_new_bool("targetbest", "if multiple top-level features "
                             "(i.e., features without parents) with exactly "
                             "one target attribute have the same target_id, "
                             "keep only the feature with the best score. If "
                             "-"TARGETGT_STRAND_OPT" is used at the same time, "
                             "this option is applied after "
                             "-"TARGETGT_STRAND_OPT".\n"
                             "Memory consumption is proportional to the input "
                             "file size(s).", &arguments->targetbest, false);
  gt_option_parser_add_option(op, option);

  /* -hascds */
  option = gt_option_new_bool("hascds", "select all top-level features which "
                              "do have a CDS child", &arguments->has_CDS,
                              false);
  gt_option_parser_add_option(op, option);

  /* -maxgenelength */
  option = gt_option_new_ulong_min("maxgenelength", "select genes up to the "
                                   "given maximum length",
                                   &arguments->max_gene_length, GT_UNDEF_ULONG,
                                   1);
  gt_option_parser_add_option(op, option);

  /* -maxgenenum */
  option = gt_option_new_ulong("maxgenenum", "select the first genes up to the "
                               "given maximum number", &arguments->max_gene_num,
                               GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, option);

  /* -mingenescore */
  option = gt_option_new_double("mingenescore", "select genes with the given "
                                "minimum score", &arguments->min_gene_score,
                                GT_UNDEF_DOUBLE);
  gt_option_parser_add_option(op, option);

  /* -maxgenescore */
  option = gt_option_new_double("maxgenescore", "select genes with the given "
                                "maximum score", &arguments->max_gene_score,
                                GT_UNDEF_DOUBLE);
  gt_option_parser_add_option(op, option);

  /* -minaveragessp */
  minaveragessp_option =
    gt_option_new_probability("minaveragessp",
                              "set the minimum average splice site probability",
                              &arguments->min_average_splice_site_prob,
                              GT_UNDEF_DOUBLE);
  gt_option_parser_add_option(op, minaveragessp_option);

  /* -singleintronfactor */
  singleintron_option =
    gt_option_new_double_min("singleintronfactor",
                             "factor to multiplicate the average splice site "
                             "probability with for single introns before "
                             "comparing it to the minimum average splice site "
                             "probability", &arguments->single_intron_factor,
                             1.0, 1.0);
  gt_option_is_development_option(singleintron_option);
  gt_option_parser_add_option(op, singleintron_option);

  /* -featurenum */
  option = gt_option_new_ulong_min("featurenum",
                                   "select feature tree occurring "
                                   "at given position in input",
                                   &arguments->feature_num, GT_UNDEF_ULONG, 1);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -filter_files */
  option = gt_option_new_filename_array("rule_files",
                                        "specify Lua files to be used "
                                        "for selection",
                                        arguments->filter_files);
  gt_option_parser_add_option(op, option);

  /* -filter_logic */
  option = gt_option_new_choice("rule_logic", "select how multiple Lua "
                                "files should be combined\nchoose from AND|OR",
                                arguments->filter_logic, filter_logic[0],
                                filter_logic);
  gt_option_parser_add_option(op, option);

  /* -nh_file */
  optiondroppedfile = gt_option_new_filename("dropped_file",
                                             "save non-selected features to "
                                             "file",
                                             arguments->dropped_file);
  gt_option_parser_add_option(op, optiondroppedfile);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* option exclusions */
  gt_option_exclude(contain_option, overlap_option);

  /* option implications */
  gt_option_imply(singleintron_option, minaveragessp_option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);

  return op;
}

static int process_gt_strand_arg(GtStr *gt_strand_char, GtStrand *strand,
                                 const char *optstr, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  if (gt_str_length(gt_strand_char)) {
    GtStrand tmpstrand = gt_strand_get(gt_str_get(gt_strand_char)[0]);
    if ((gt_str_length(gt_strand_char) > 1) ||
        (tmpstrand == GT_NUM_OF_STRAND_TYPES)) {
      gt_error_set(err, "argument to option -%s must be one of '"
                        GT_STRAND_CHARS"'", optstr);
      had_err = -1;
    }
    if (!had_err)
      *strand = tmpstrand;
  }
  return had_err;
}

static int gt_select_arguments_check(GT_UNUSED int rest_argc,
                                     void *tool_arguments, GtError *err)
{
  SelectArguments *arguments = tool_arguments;
  int had_err;
  gt_error_check(err);
  gt_assert(arguments);
  had_err = process_gt_strand_arg(arguments->gt_strand_char, &arguments->strand,
                                  GT_STRAND_OPT, err);
  if (!had_err) {
    had_err = process_gt_strand_arg(arguments->targetgt_strand_char,
                                    &arguments->targetstrand,
                                    TARGETGT_STRAND_OPT, err);
  }

  return had_err;
}

static int print_to_file_drophandler(GtGenomeNode *gn, void *data,
                                     GtError *err)
{
  GtNodeVisitor *v;
  gt_assert(gn && data);
  v = (GtNodeVisitor*) data;
  return gt_genome_node_accept(gn, v, err);
}

static int default_drophandler(GT_UNUSED GtGenomeNode *gn, GT_UNUSED void *data,
                               GT_UNUSED GtError *err)
{
  gt_assert(gn);
  return 0;
}

static int gt_select_runner(int argc, const char **argv, int parsed_args,
                            void *tool_arguments, GtError *err)
{
  SelectArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream, *select_stream,
               *targetbest_select_stream = NULL, *gff3_out_stream;
  int had_err;
  GtFile *drop_file = NULL;
  GtNodeVisitor *gff3outvis = NULL;
  gt_error_check(err);
  gt_assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  if (arguments->verbose && arguments->outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create a filter stream */
  select_stream = gt_select_stream_new(gff3_in_stream, arguments->seqid,
                                       arguments->source,
                                       &arguments->contain_range,
                                       &arguments->overlap_range,
                                       arguments->strand,
                                       arguments->targetstrand,
                                       arguments->has_CDS,
                                       arguments->max_gene_length,
                                       arguments->max_gene_num,
                                       arguments->min_gene_score,
                                       arguments->max_gene_score,
                                       arguments->min_average_splice_site_prob,
                                       arguments->feature_num,
                                       arguments->filter_files,
                                       arguments->filter_logic,
                                       err);

  if (select_stream) {
    GtSelectStream *fs = (GtSelectStream*) select_stream;

    if (gt_str_length(arguments->dropped_file) > 0) {
      drop_file = gt_file_new(gt_str_get(arguments->dropped_file), "w", err);
      gff3outvis = gt_gff3_visitor_new(drop_file);
      gt_select_stream_set_drophandler(fs, print_to_file_drophandler,
                                       (void*) gff3outvis);
    } else {
      gt_select_stream_set_drophandler(fs, default_drophandler, NULL);
    }

    gt_select_stream_set_single_intron_factor(select_stream,
                                              arguments->single_intron_factor);

    if (arguments->targetbest)
      targetbest_select_stream = gt_targetbest_select_stream_new(select_stream);

    /* create a gff3 output stream */
    gff3_out_stream = gt_gff3_out_stream_new(arguments->targetbest
                                             ? targetbest_select_stream
                                             : select_stream,
                                             arguments->outfp);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(gff3_out_stream, err);

    /* free */
    gt_node_stream_delete(gff3_out_stream);
    gt_node_stream_delete(select_stream);
    gt_node_stream_delete(targetbest_select_stream);
  } else {
    had_err = -1;
  }
  gt_file_delete(drop_file);
  gt_node_visitor_delete(gff3outvis);
  gt_node_stream_delete(gff3_in_stream);
  return had_err;
}

GtTool* gt_select(void)
{
  return gt_tool_new(gt_select_arguments_new,
                     gt_select_arguments_delete,
                     gt_select_option_parser_new,
                     gt_select_arguments_check,
                     gt_select_runner);
}
