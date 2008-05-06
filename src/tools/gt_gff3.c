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

#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/outputfile.h"
#include "libgtcore/undef.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/add_introns_stream.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/mergefeat_stream_sorted.h"
#include "libgtext/sort_stream.h"
#include "tools/gt_gff3.h"

typedef struct {
  bool sort,
       checkids,
       mergefeat,
       addintrons,
       verbose;
  long offset;
  Str *offsetfile;
  OutputFileInfo *ofi;
  GenFile *outfp;
} GFF3Arguments;

static void* gt_gff3_arguments_new(void)
{
  GFF3Arguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->offsetfile = str_new();
  arguments->ofi = outputfileinfo_new();
  return arguments;
}

static void gt_gff3_arguments_delete(void *tool_arguments)
{
  GFF3Arguments *arguments = tool_arguments;
  if (!arguments) return;
  genfile_close(arguments->outfp);
  outputfileinfo_delete(arguments->ofi);
  str_delete(arguments->offsetfile);
  ma_free(arguments);
}

static OptionParser* gt_gff3_option_parser_new(void *tool_arguments)
{
  GFF3Arguments *arguments = tool_arguments;
  OptionParser *op;
  Option *sort_option, *mergefeat_option, *addintrons_option, *offset_option,
         *offsetfile_option, *option;
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

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  /* set comment function */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);

  return op;
}

static int gt_gff3_runner(int argc, const char **argv, int parsed_args,
                          void *tool_arguments, Error *err)
{
  GenomeStream *gff3_in_stream,
               *sort_stream = NULL,
               *mergefeat_stream = NULL,
               *add_introns_stream = NULL,
               *gff3_out_stream = NULL,
               *last_stream;
  GFF3Arguments *arguments = tool_arguments;
  GenomeNode *gn;
  int had_err = 0;

  error_check(err);
  assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments->verbose &&
                                               arguments->outfp,
                                               arguments->checkids);
  last_stream = gff3_in_stream;

  /* set offset (if necessary) */
  if (arguments->offset != UNDEF_LONG)
    gff3_in_stream_set_offset(gff3_in_stream, arguments->offset);

  /* set offsetfile (if necessary) */
  if (str_length(arguments->offsetfile)) {
    had_err = gff3_in_stream_set_offsetfile(gff3_in_stream,
                                            arguments->offsetfile, err);
  }

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
  if (!had_err)
    gff3_out_stream = gff3_out_stream_new(last_stream, arguments->outfp);

  /* pull the features through the stream and free them afterwards */
  if (!had_err) {
    while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, err)) &&
           gn) {
      genome_node_rec_delete(gn);
    }
  }

  /* free */
  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(sort_stream);
  genome_stream_delete(mergefeat_stream);
  genome_stream_delete(add_introns_stream);
  genome_stream_delete(gff3_in_stream);

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
