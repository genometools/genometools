/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/mergefeat_stream_sorted.h"
#include "extended/sort_stream.h"
#include "extended/type_factory_builtin.h"
#include "extended/type_factory_obo.h"
#include "tools/gt_gff3.h"

typedef struct {
  bool sort,
       checkids,
       mergefeat,
       addintrons,
       verbose,
       typecheck_built_in,
       tidy;
  long offset;
  GT_Str *offsetfile,
      *typecheck;
  unsigned long width;
  OutputFileInfo *ofi;
  GT_GenFile *outfp;
} GFF3Arguments;

static void* gt_gff3_arguments_new(void)
{
  GFF3Arguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->offsetfile = gt_str_new();
  arguments->typecheck = gt_str_new();
  arguments->ofi = outputfileinfo_new();
  return arguments;
}

static void gt_gff3_arguments_delete(void *tool_arguments)
{
  GFF3Arguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_genfile_close(arguments->outfp);
  outputfileinfo_delete(arguments->ofi);
  gt_str_delete(arguments->typecheck);
  gt_str_delete(arguments->offsetfile);
  gt_free(arguments);
}

static OptionParser* gt_gff3_option_parser_new(void *tool_arguments)
{
  GFF3Arguments *arguments = tool_arguments;
  OptionParser *op;
  Option *sort_option, *mergefeat_option, *addintrons_option, *offset_option,
         *offsetfile_option, *typecheck_option, *built_in_option, *option;
  assert(arguments);

  /* init */
  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Parse, possibly transform, and output GFF3 files.");

  /* -sort */
  sort_option = option_new_bool("sort", "sort the GFF3 features (memory "
                                "consumption is O(file_size))",
                                &arguments->sort, false);
  option_parser_add_option(op, sort_option);

  /* -checkids */
  option = option_new_bool("checkids", "make sure the ID attributes are unique "
                           "within the scope of each GFF3_file, as required by "
                           "GFF3 specification\n"
                           "(memory consumption is O(file_size))",
                           &arguments->checkids, false);
  option_parser_add_option(op, option);

  /* -mergefeat */
  mergefeat_option = option_new_bool("mergefeat", "merge adjacent features of "
                                     "the same type", &arguments->mergefeat,
                                     false);
  option_is_development_option(mergefeat_option);
  option_imply(mergefeat_option, sort_option);
  option_parser_add_option(op, mergefeat_option);

  /* -addintrons */
  addintrons_option = option_new_bool("addintrons", "add intron features "
                                      "between existing exon features",
                                      &arguments->addintrons, false);
  option_parser_add_option(op, addintrons_option);

  /* -offset */
  offset_option = option_new_long("offset",
                                 "transform all features by the given offset",
                                  &arguments->offset, UNDEF_LONG);
  option_parser_add_option(op, offset_option);

  /* -offsetfile */
  offsetfile_option = option_new_filename("offsetfile", "transform all "
                                          "features by the offsets given in "
                                          "file", arguments->offsetfile);
  option_parser_add_option(op, offsetfile_option);
  option_exclude(offset_option, offsetfile_option);

  /* -typecheck */
  typecheck_option = option_new_filename("typecheck", "check GFF3 types "
                                         "against \"id\" and \"name\" tags "
                                         "in given OBO file",
                                         arguments->typecheck);
  option_parser_add_option(op, typecheck_option);

  /* -typecheck-built-in */
  built_in_option = option_new_bool("typecheck-built-in", "use built-in type "
                                    "checker", &arguments->typecheck_built_in,
                                    false);
  option_is_development_option(built_in_option);
  option_parser_add_option(op, built_in_option);
  option_exclude(typecheck_option, built_in_option);

  /* -tidy */
  option = option_new_bool("tidy", "try to tidy the file", &arguments->tidy,
                           false);
  option_is_development_option(option);
  option_parser_add_option(op, option);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* -width */
  option = option_new_ulong("width", "set output width for showing of embedded "
                            "FASTA sequences\n(0 disables formatting)",
                            &arguments->width, 0);
  option_parser_add_option(op, option);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  /* set comment function */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);

  return op;
}

static int gt_gff3_runner(int argc, const char **argv, int parsed_args,
                          void *tool_arguments, GT_Error *err)
{
  GT_TypeFactory *ftf = NULL;
  GenomeStream *gff3_in_stream,
               *sort_stream = NULL,
               *mergefeat_stream = NULL,
               *add_introns_stream = NULL,
               *gff3_out_stream = NULL,
               *last_stream;
  GFF3Arguments *arguments = tool_arguments;
  GT_GenomeNode *gn;
  int had_err = 0;

  gt_error_check(err);
  assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments->verbose &&
                                               arguments->outfp,
                                               arguments->checkids);
  last_stream = gff3_in_stream;

  /* set different type checker if necessary */
  if (arguments->typecheck_built_in) {
      ftf = gt_type_factory_builtin_new();
      gff3_in_stream_set_feature_type_factory(gff3_in_stream, ftf);
  }
  if (gt_str_length(arguments->typecheck)) {
    if (!(ftf = gt_type_factory_obo_new(gt_str_get(arguments->typecheck),
                                             err))) {
        had_err = -1;
    }
    if (!had_err)
      gff3_in_stream_set_feature_type_factory(gff3_in_stream, ftf);
  }

  /* set offset (if necessary) */
  if (!had_err && arguments->offset != UNDEF_LONG)
    gff3_in_stream_set_offset(gff3_in_stream, arguments->offset);

  /* set offsetfile (if necessary) */
  if (!had_err && gt_str_length(arguments->offsetfile)) {
    had_err = gff3_in_stream_set_offsetfile(gff3_in_stream,
                                            arguments->offsetfile, err);
  }

  /* enable tidy mode (if necessary) */
  if (!had_err && arguments->tidy)
    gff3_in_stream_enable_tidy_mode(gff3_in_stream);

  /* create sort stream (if necessary) */
  if (!had_err && arguments->sort) {
    sort_stream = sort_stream_new(gff3_in_stream);
    last_stream = sort_stream;
  }

  /* create merge feature stream (if necessary) */
  if (!had_err && arguments->mergefeat) {
    assert(sort_stream);
    mergefeat_stream = mergefeat_stream_sorted_new(sort_stream);
    last_stream = mergefeat_stream;
  }

  /* create addintrons stream (if necessary) */
  if (!had_err && arguments->addintrons) {
    assert(last_stream);
    add_introns_stream = add_introns_stream_new(last_stream);
    last_stream = add_introns_stream;
  }

  /* create gff3 output stream */
  if (!had_err) {
    gff3_out_stream = gff3_out_stream_new(last_stream, arguments->outfp);
    gff3_out_stream_set_fasta_width(gff3_out_stream, arguments->width);
  }

  /* pull the features through the stream and free them afterwards */
  if (!had_err) {
    while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, err)) &&
           gn) {
      gt_genome_node_rec_delete(gn);
    }
  }

  /* free */
  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(sort_stream);
  genome_stream_delete(mergefeat_stream);
  genome_stream_delete(add_introns_stream);
  genome_stream_delete(gff3_in_stream);
  gt_type_factory_delete(ftf);

  return had_err;
}

Tool* gt_gff3(void)
{
  return tool_new(gt_gff3_arguments_new,
                  gt_gff3_arguments_delete,
                  gt_gff3_option_parser_new,
                  NULL,
                  gt_gff3_runner);
}
