/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/ma.h"
#include "core/option.h"
#include "core/outputfile.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/filter_stream.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream.h"
#include "extended/targetbest_filter_stream.h"
#include "tools/gt_filter.h"

#define GT_STRAND_OPT  "strand"
#define TARGETGT_STRAND_OPT  "targetstrand"

typedef struct {
  bool verbose,
       has_CDS,
       targetbest;
  GtStr *seqid,
        *typefilter,
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
         min_average_splice_site_prob;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} FilterArguments;

static void* gt_filter_arguments_new(void)
{
  FilterArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->seqid = gt_str_new();
  arguments->typefilter = gt_str_new();
  arguments->gt_strand_char = gt_str_new();
  arguments->strand = GT_NUM_OF_STRAND_TYPES;
  arguments->targetgt_strand_char = gt_str_new();
  arguments->targetstrand = GT_NUM_OF_STRAND_TYPES;
  arguments->ofi = gt_outputfileinfo_new();
  return arguments;
}

static void gt_filter_arguments_delete(void *tool_arguments)
{
  FilterArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_str_delete(arguments->targetgt_strand_char);
  gt_str_delete(arguments->gt_strand_char);
  gt_str_delete(arguments->typefilter);
  gt_str_delete(arguments->seqid);
  gt_free(arguments);
}

static GtOptionParser* gt_filter_option_parser_new(void *tool_arguments)
{
  FilterArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *contain_option, *overlap_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                            "Filter GFF3 files.");

  /* -seqid */
  option = gt_option_new_string("seqid", "seqid a feature must have to pass "
                                "the filter (excluding comments)",
                                arguments->seqid,
                                NULL);
  gt_option_parser_add_option(op, option);

  /* -typefilter */
  option = gt_option_new_string("typefilter", "filter out all features of the "
                                "given type", arguments->typefilter, NULL);
  /* XXX */
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -contain */
  contain_option = gt_option_new_range("contain", "filter out all features "
                                       "which are not contained in the given "
                                       "range",
                                       &arguments->contain_range, NULL);
  gt_option_parser_add_option(op, contain_option);

  /* -overlap */
  overlap_option = gt_option_new_range("overlap", "filter out all features "
                                       "which do not overlap with the given "
                                       "range",
                                       &arguments->overlap_range, NULL);
  gt_option_parser_add_option(op, overlap_option);

  /* -strand */
  option = gt_option_new_string(GT_STRAND_OPT, "filter out all top-level "
                                "features (i.e., features without parents) "
                                "whose strand is different from the given one "
                                "(must be one of '"
                                GT_STRAND_CHARS"')", arguments->gt_strand_char,
                                NULL);
  gt_option_parser_add_option(op, option);

  /* -targetstrand */
  option = gt_option_new_string(TARGETGT_STRAND_OPT, "filter out all top-level "
                             "features (i.e., features without parents) which "
                             "have exactly one target attribute whose strand "
                             "is different from the given one (must be one of '"
                             GT_STRAND_CHARS"')",
                             arguments->targetgt_strand_char, NULL);
  gt_option_parser_add_option(op, option);

  /* -targetbest */
  option = gt_option_new_bool("targetbest", "if multiple top-level features "
                           "(i.e., features without parents) with exactly one "
                           "target attribute have the same target_id, keep "
                           "only the feature with the best score. If "
                           "-"TARGETGT_STRAND_OPT" is used at the same time, "
                           "this option is applied after "
                           "-"TARGETGT_STRAND_OPT".\n"
                           "Memory consumption is O(file_size).",
                           &arguments->targetbest, false);
  gt_option_parser_add_option(op, option);

  /* -hascds */
  option = gt_option_new_bool("hascds", "filter out all top-level features "
                              "which do not have a CDS child",
                              &arguments->has_CDS,
                              false);
  gt_option_parser_add_option(op, option);

  /* -maxgenelength */
  option = gt_option_new_ulong_min("maxgenelength", "the maximum length a gene "
                                "can have to pass the filter",
                                &arguments->max_gene_length, UNDEF_ULONG, 1);
  gt_option_parser_add_option(op, option);

  /* -maxgenenum */
  option = gt_option_new_ulong("maxgenenum", "the maximum number of genes "
                               "which can pass the filter",
                               &arguments->max_gene_num,
                               UNDEF_ULONG);
  gt_option_parser_add_option(op, option);

  /* -mingenescore */
  option = gt_option_new_double("mingenescore", "the minimum score a gene must "
                                "have to pass the filter",
                                &arguments->min_gene_score, UNDEF_DOUBLE);
  gt_option_parser_add_option(op, option);

  /* -maxgenescore */
  option = gt_option_new_double("maxgenescore", "the maximum score a gene can "
                                "have to pass the filter",
                                &arguments->max_gene_score, UNDEF_DOUBLE);
  gt_option_parser_add_option(op, option);

  /* -minaveragessp */
  option = gt_option_new_probability("minaveragessp", "set the minimum average "
                                     "splice site probability",
                                     &arguments->min_average_splice_site_prob,
                                     UNDEF_DOUBLE);
  gt_option_parser_add_option(op, option);

  /* -featurenum */
  option = gt_option_new_ulong_min("featurenum",
                                   "select feature tree occurring "
                                   "at given position in input",
                                   &arguments->feature_num, UNDEF_ULONG, 1);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* option exclusions */
  gt_option_exclude(contain_option, overlap_option);

  /* output file options */
  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);

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

static int gt_filter_arguments_check(GT_UNUSED int rest_argc,
                                     void *tool_arguments, GtError *err)
{
  FilterArguments *arguments = tool_arguments;
  int had_err;
  gt_error_check(err);
  gt_assert(arguments);
  had_err = process_gt_strand_arg(arguments->gt_strand_char, &arguments->strand,
                               GT_STRAND_OPT, err);
  if (!had_err) {
    had_err = process_gt_strand_arg(arguments->targetgt_strand_char,
                                 &arguments->targetstrand, TARGETGT_STRAND_OPT,
                                 err);
  }
  return had_err;
}

static int gt_filter_runner(int argc, const char **argv, int parsed_args,
                           void *tool_arguments, GtError *err)
{
  FilterArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream, *filter_stream,
               *targetbest_filter_stream = NULL, *gff3_out_stream;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  if (arguments->verbose && arguments->outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create a filter stream */
  filter_stream = gt_filter_stream_new(gff3_in_stream, arguments->seqid,
                                       arguments->typefilter,
                                       arguments->contain_range,
                                       arguments->overlap_range,
                                       arguments->strand,
                                       arguments->targetstrand,
                                       arguments->has_CDS,
                                       arguments->max_gene_length,
                                       arguments->max_gene_num,
                                       arguments->min_gene_score,
                                       arguments->max_gene_score,
                                       arguments->min_average_splice_site_prob,
                                       arguments->feature_num);

  if (arguments->targetbest)
    targetbest_filter_stream = gt_targetbest_filter_stream_new(filter_stream);

  /* create a gff3 output stream */
  gff3_out_stream = gt_gff3_out_stream_new(arguments->targetbest
                                           ? targetbest_filter_stream
                                           : filter_stream,
                                           arguments->outfp);

  /* pull the features through the stream and free them afterwards */
  had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(filter_stream);
  gt_node_stream_delete(targetbest_filter_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_filter(void)
{
  return gt_tool_new(gt_filter_arguments_new,
                  gt_filter_arguments_delete,
                  gt_filter_option_parser_new,
                  gt_filter_arguments_check,
                  gt_filter_runner);
}
