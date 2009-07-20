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

#include <string.h>
#include "core/ma.h"
#include "core/option.h"
#include "core/outputfile.h"
#include "core/undef.h"
#include "core/versionfunc.h"
#include "extended/add_introns_stream.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/merge_feature_stream.h"
#include "extended/sort_stream.h"
#include "extended/type_checker_builtin.h"
#include "extended/type_checker_obo.h"
#include "tools/gt_gff3.h"

typedef struct {
  bool sort,
       checkids,
       retainids,
       mergefeat,
       addintrons,
       verbose,
       typecheck_built_in,
       tidy;
  long offset;
  GtStr *offsetfile,
        *typecheck;
  unsigned long width;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GFF3Arguments;

static void* gt_gff3_arguments_new(void)
{
  GFF3Arguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->offsetfile = gt_str_new();
  arguments->typecheck = gt_str_new();
  arguments->ofi = gt_outputfileinfo_new();
  return arguments;
}

static void gt_gff3_arguments_delete(void *tool_arguments)
{
  GFF3Arguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_str_delete(arguments->typecheck);
  gt_str_delete(arguments->offsetfile);
  gt_free(arguments);
}

static GtOptionParser* gt_gff3_option_parser_new(void *tool_arguments)
{
  GFF3Arguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *sort_option, *mergefeat_option, *addintrons_option, *offset_option,
         *offsetfile_option, *typecheck_option, *built_in_option, *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                         "Parse, possibly transform, and output GFF3 files.");

  /* -sort */
  sort_option = gt_option_new_bool("sort", "sort the GFF3 features (memory "
                                "consumption is O(file_size))",
                                &arguments->sort, false);
  gt_option_parser_add_option(op, sort_option);

  /* -tidy */
  option = gt_option_new_bool("tidy", "try to tidy the GFF3 files up during "
                           "parsing", &arguments->tidy, false);
  gt_option_parser_add_option(op, option);

  /* -retainids */
  option = gt_option_new_bool("retainids",
                           "when available, use the original IDs provided"
                           "in the source file\n"
                           "(memory consumption is O(file_size))",
                           &arguments->retainids, false);
  gt_option_parser_add_option(op, option);

  /* -checkids */
  option = gt_option_new_bool("checkids",
                           "make sure the ID attributes are unique "
                           "within the scope of each GFF3_file, as required by "
                           "GFF3 specification\n"
                           "(memory consumption is O(file_size))",
                           &arguments->checkids, false);
  gt_option_parser_add_option(op, option);

  /* -mergefeat */
  mergefeat_option = gt_option_new_bool("mergefeat",
                                     "merge adjacent features of "
                                     "the same type", &arguments->mergefeat,
                                     false);
  gt_option_is_development_option(mergefeat_option);
  gt_option_imply(mergefeat_option, sort_option);
  gt_option_parser_add_option(op, mergefeat_option);

  /* -addintrons */
  addintrons_option = gt_option_new_bool("addintrons", "add intron features "
                                      "between existing exon features",
                                      &arguments->addintrons, false);
  gt_option_parser_add_option(op, addintrons_option);

  /* -offset */
  offset_option = gt_option_new_long("offset",
                                 "transform all features by the given offset",
                                  &arguments->offset, UNDEF_LONG);
  gt_option_parser_add_option(op, offset_option);

  /* -offsetfile */
  offsetfile_option = gt_option_new_filename("offsetfile", "transform all "
                                          "features by the offsets given in "
                                          "file", arguments->offsetfile);
  gt_option_parser_add_option(op, offsetfile_option);
  gt_option_exclude(offset_option, offsetfile_option);

  /* -typecheck */
  typecheck_option = gt_option_new_filename("typecheck", "check GFF3 types "
                                         "against \"id\" and \"name\" tags "
                                         "in given OBO file",
                                         arguments->typecheck);
  gt_option_parser_add_option(op, typecheck_option);

  /* -typecheck-built-in */
  built_in_option = gt_option_new_bool("typecheck-built-in",
                                    "use built-in type "
                                    "checker", &arguments->typecheck_built_in,
                                    false);
  gt_option_is_development_option(built_in_option);
  gt_option_parser_add_option(op, built_in_option);
  gt_option_exclude(typecheck_option, built_in_option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* -width */
  option = gt_option_new_ulong("width",
                            "set output width for showing of embedded "
                            "FASTA sequences\n(0 disables formatting)",
                            &arguments->width, 0);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  /* set comment function */
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);

  return op;
}

static int gt_gff3_runner(int argc, const char **argv, int parsed_args,
                          void *tool_arguments, GtError *err)
{
  GtTypeChecker *type_checker = NULL;
  GtNodeStream *gff3_in_stream,
               *sort_stream = NULL,
               *merge_feature_stream = NULL,
               *add_introns_stream = NULL,
               *gff3_out_stream = NULL,
               *last_stream;
  GFF3Arguments *arguments = tool_arguments;
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

  last_stream = gff3_in_stream;

  /* set different type checker if necessary */
  if (arguments->typecheck_built_in) {
      type_checker = gt_type_checker_builtin_new();
      gt_gff3_in_stream_set_type_checker(gff3_in_stream, type_checker);
  }
  if (gt_str_length(arguments->typecheck)) {
    type_checker = gt_type_checker_obo_new(gt_str_get(arguments->typecheck),
                                           err);
    if (!type_checker)
      had_err = -1;
    if (!had_err)
      gt_gff3_in_stream_set_type_checker(gff3_in_stream, type_checker);
  }

  /* set offset (if necessary) */
  if (!had_err && arguments->offset != UNDEF_LONG)
    gt_gff3_in_stream_set_offset(gff3_in_stream, arguments->offset);

  /* set offsetfile (if necessary) */
  if (!had_err && gt_str_length(arguments->offsetfile)) {
    had_err = gt_gff3_in_stream_set_offsetfile(gff3_in_stream,
                                               arguments->offsetfile, err);
  }

  /* enable tidy mode (if necessary) */
  if (!had_err && arguments->tidy)
    gt_gff3_in_stream_enable_tidy_mode(gff3_in_stream);

  /* create sort stream (if necessary) */
  if (!had_err && arguments->sort) {
    sort_stream = gt_sort_stream_new(gff3_in_stream);
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

  /* create gff3 output stream */
  if (!had_err) {
    gff3_out_stream = gt_gff3_out_stream_new(last_stream, arguments->outfp);
    gt_gff3_out_stream_set_fasta_width(gff3_out_stream, arguments->width);
  }

  if (!had_err && arguments->retainids)
    gt_gff3_out_stream_retain_id_attributes(gff3_out_stream);

  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(sort_stream);
  gt_node_stream_delete(merge_feature_stream);
  gt_node_stream_delete(add_introns_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_type_checker_delete(type_checker);

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
