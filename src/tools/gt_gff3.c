/*
  Copyright (c) 2005-2013 Gordon Gremme <gordon@gremme.org>
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
#include "core/versionfunc.h"
#include "extended/add_introns_stream_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gff3_linesorted_out_stream.h"
#include "extended/gff3_parser.h"
#include "extended/gtdatahelp.h"
#include "extended/load_stream.h"
#include "extended/merge_feature_stream_api.h"
#include "extended/set_source_visitor_api.h"
#include "extended/sort_stream_api.h"
#include "extended/typecheck_info.h"
#include "extended/visitor_stream_api.h"
#include "extended/xrfcheck_info.h"
#include "tools/gt_gff3.h"

typedef struct {
  bool sort,
       sortlines,
       load,
       retainids,
       checkids,
       addids,
       mergefeat,
       addintrons,
       verbose,
       strict,
       tidy,
       show,
       fixboundaries;
  GtWord offset;
  GtStr *offsetfile, *newsource;
  GtUword width;
  GtTypecheckInfo *tci;
  GtXRFCheckInfo *xci;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GFF3Arguments;

static void* gt_gff3_arguments_new(void)
{
  GFF3Arguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->newsource = gt_str_new();
  arguments->offsetfile = gt_str_new();
  arguments->tci = gt_typecheck_info_new();
  arguments->xci = gt_xrfcheck_info_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_gff3_arguments_delete(void *tool_arguments)
{
  GFF3Arguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_str_delete(arguments->newsource);
  gt_output_file_info_delete(arguments->ofi);
  gt_typecheck_info_delete(arguments->tci);
  gt_xrfcheck_info_delete(arguments->xci);
  gt_str_delete(arguments->offsetfile);
  gt_free(arguments);
}

static GtOptionParser* gt_gff3_option_parser_new(void *tool_arguments)
{
  GFF3Arguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *sort_option, *load_option, *strict_option, *tidy_option,
           *mergefeat_option, *addintrons_option, *offset_option,
           *offsetfile_option, *setsource_option, *sortlines_option, *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file ...]", "Parse, possibly "
                            "transform, and output GFF3 files.");

  /* -sort */
  sort_option = gt_option_new_bool("sort", "sort the GFF3 features (memory "
                                   "consumption is proportional to the input "
                                   "file size(s))",
                                   &arguments->sort, false);
  gt_option_parser_add_option(op, sort_option);

  /* -sortlines */
  sortlines_option = gt_option_new_bool("sortlines", "sort the GFF3 features "
                                        "on a strict line basis (not sorted as"
                                        "defined by GenomeTools)",
                                        &arguments->sortlines, false);
  gt_option_imply(sortlines_option, sort_option);
  gt_option_parser_add_option(op, sortlines_option);

  /* -strict */
  strict_option = gt_option_new_bool("strict", "be very strict during GFF3 "
                                     "parsing (stricter than the specification "
                                     "requires)", &arguments->strict, false);
  gt_option_is_development_option(strict_option);
  gt_option_parser_add_option(op, strict_option);

  /* -tidy */
  tidy_option = gt_option_new_bool("tidy", "try to tidy the GFF3 files up "
                                   "during parsing", &arguments->tidy, false);
  gt_option_parser_add_option(op, tidy_option);
  gt_option_exclude(strict_option, tidy_option);

  /* -retainids */
  option = gt_option_new_bool("retainids",
                              "when available, use the original IDs provided "
                              "in the source file\n"
                              "(memory consumption is proportional to the "
                              "input file size(s))", &arguments->retainids,
                              false);
  gt_option_parser_add_option(op, option);

  /* -checkids */
  option = gt_option_new_bool("checkids",
                              "make sure the ID attributes are unique "
                              "within the scope of each GFF3_file, as required "
                              "by GFF3 specification\n"
                              "(memory consumption is proportional to the "
                              "input file size(s))", &arguments->checkids,
                              false);
  gt_option_parser_add_option(op, option);

  /* -addids */
  option = gt_option_new_bool("addids", "add missing \""
                              GT_GFF_SEQUENCE_REGION"\" lines automatically",
                              &arguments->addids, true);
  gt_option_parser_add_option(op, option);

  /* -fixregionboundaries */
  option = gt_option_new_bool("fixregionboundaries", "automatically adjust \""
                              GT_GFF_SEQUENCE_REGION"\" lines to contain all "
                              "their features (memory consumption is "
                              "proportional to the input file size(s))",
                              &arguments->fixboundaries, false);
  gt_option_parser_add_option(op, option);

  /* -mergefeat */
  mergefeat_option = gt_option_new_bool("mergefeat",
                                        "merge adjacent features of the same "
                                        "type", &arguments->mergefeat, false);
  gt_option_is_development_option(mergefeat_option);
  gt_option_imply(mergefeat_option, sort_option);
  gt_option_parser_add_option(op, mergefeat_option);

  /* -load */
  load_option = gt_option_new_bool("load", "load the GFF3 features into memory "
                                   "(requires space proportional to the input "
                                   "file size(s))",
                                   &arguments->load, false);
  gt_option_is_development_option(load_option);
  gt_option_parser_add_option(op, load_option);

  /* -addintrons */
  addintrons_option = gt_option_new_bool("addintrons", "add intron features "
                                         "between existing exon features",
                                         &arguments->addintrons, false);
  gt_option_parser_add_option(op, addintrons_option);

  /* -offset */
  offset_option = gt_option_new_word("offset", "transform all features by the "
                                     "given offset", &arguments->offset,
                                     GT_UNDEF_WORD);
  gt_option_parser_add_option(op, offset_option);

  /* -offsetfile */
  offsetfile_option = gt_option_new_filename("offsetfile", "transform all "
                                             "features by the offsets given in "
                                             "file", arguments->offsetfile);
  gt_option_parser_add_option(op, offsetfile_option);
  gt_option_exclude(offset_option, offsetfile_option);

  /* -setsource */
  setsource_option = gt_option_new_string("setsource", "set the 'source' "
                                          "value (2nd column) of each feature",
                                          arguments->newsource, NULL);
  gt_option_parser_add_option(op, setsource_option);

  /* typecheck options */
  gt_typecheck_info_register_options(arguments->tci, op);

  /* xrfcheck options */
  gt_xrfcheck_info_register_options(arguments->xci, op);

  /* -show */
  option = gt_option_new_bool("show", "show GFF3 output", &arguments->show,
                              true);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* -width */
  option = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  /* set comment function */
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);

  return op;
}

static int gt_gff3_runner(int argc, const char **argv, int parsed_args,
                          void *tool_arguments, GtError *err)
{
  GFF3Arguments *arguments = tool_arguments;
  GtTypeChecker *type_checker = NULL;
  GtXRFChecker *xrf_checker = NULL;
  GtNodeStream *gff3_in_stream,
               *sort_stream = NULL,
               *load_stream = NULL,
               *merge_feature_stream = NULL,
               *add_introns_stream = NULL,
               *set_source_stream = NULL,
               *gff3_out_stream = NULL,
               *last_stream;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  if (arguments->verbose && arguments->outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);
  if (arguments->checkids)
    gt_gff3_in_stream_check_id_attributes((GtGFF3InStream*) gff3_in_stream);
  if (!arguments->addids)
    gt_gff3_in_stream_disable_add_ids(gff3_in_stream);

  last_stream = gff3_in_stream;

  /* set different type checker if necessary */
  if (gt_typecheck_info_option_used(arguments->tci)) {
    type_checker = gt_typecheck_info_create_type_checker(arguments->tci, err);
    if (!type_checker)
      had_err = -1;
    if (!had_err)
      gt_gff3_in_stream_set_type_checker(gff3_in_stream, type_checker);
  }

  /* set XRF checker if necessary */
  if (gt_xrfcheck_info_option_used(arguments->xci)) {
    xrf_checker = gt_xrfcheck_info_create_xrf_checker(arguments->xci, err);
    if (!xrf_checker)
      had_err = -1;
    if (!had_err)
      gt_gff3_in_stream_set_xrf_checker(gff3_in_stream, xrf_checker);
  }

  /* set offset (if necessary) */
  if (!had_err && arguments->offset != GT_UNDEF_WORD)
    gt_gff3_in_stream_set_offset(gff3_in_stream, arguments->offset);

  /* set offsetfile (if necessary) */
  if (!had_err && gt_str_length(arguments->offsetfile)) {
    had_err = gt_gff3_in_stream_set_offsetfile(gff3_in_stream,
                                               arguments->offsetfile, err);
  }

  /* enable strict mode (if necessary) */
  if (!had_err && arguments->strict)
    gt_gff3_in_stream_enable_strict_mode((GtGFF3InStream*) gff3_in_stream);
  /* enable tidy mode (if necessary) */
  if (!had_err && arguments->tidy)
    gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream*) gff3_in_stream);

  if (!had_err && arguments->fixboundaries)
    gt_gff3_in_stream_fix_region_boundaries((GtGFF3InStream*) gff3_in_stream);

  /* create load stream (if necessary) */
  if (!had_err && arguments->load) {
    load_stream = gt_load_stream_new(last_stream);
    last_stream = load_stream;
  }

  /* create sort stream (if necessary) */
  if (!had_err && arguments->sort) {
    sort_stream = gt_sort_stream_new(last_stream);
    last_stream = sort_stream;
  }

  /* create merge feature stream (if necessary) */
  if (!had_err && arguments->mergefeat) {
    gt_assert(sort_stream);
    merge_feature_stream = gt_merge_feature_stream_new(sort_stream);
    last_stream = merge_feature_stream;
  }

  /* create addintrons stream (if necessary) */
  if (!had_err && arguments->addintrons) {
    gt_assert(last_stream);
    add_introns_stream = gt_add_introns_stream_new(last_stream);
    last_stream = add_introns_stream;
  }

  /* create setsource stream (if necessary) */
  if (!had_err && gt_str_length(arguments->newsource) > 0) {
    gt_assert(last_stream);
    GtNodeVisitor *ssv = gt_set_source_visitor_new(arguments->newsource);
    set_source_stream = gt_visitor_stream_new(last_stream, ssv);
    last_stream = set_source_stream;
  }

  /* create gff3 output stream */
  if (!had_err && arguments->show) {
    if (arguments->sortlines) {
      gff3_out_stream = gt_gff3_linesorted_out_stream_new(last_stream,
                                                          arguments->outfp);
    } else {
      gff3_out_stream = gt_gff3_out_stream_new(last_stream, arguments->outfp);
      gt_gff3_out_stream_set_fasta_width((GtGFF3OutStream*) gff3_out_stream,
                                         arguments->width);
      if (arguments->retainids)
        gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream*)
                                                               gff3_out_stream);
    }
    gt_assert(gff3_out_stream);
    last_stream = gff3_out_stream;
  }

  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(last_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(sort_stream);
  gt_node_stream_delete(load_stream);
  gt_node_stream_delete(merge_feature_stream);
  gt_node_stream_delete(add_introns_stream);
  gt_node_stream_delete(set_source_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_type_checker_delete(type_checker);
  gt_xrf_checker_delete(xrf_checker);

  return had_err;
}

GtTool* gt_gff3(void)
{
  return gt_tool_new(gt_gff3_arguments_new,
                     gt_gff3_arguments_delete,
                     gt_gff3_option_parser_new,
                     NULL,
                     gt_gff3_runner);
}
