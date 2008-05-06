/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/unused.h"
#include "libgtext/cds_stream.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/seqid2file.h"
#include "tools/gt_cds.h"

#define GT_CDS_SOURCE_TAG "gt cds"

typedef struct {
  bool verbose;
  Str *seqfile,
      *regionmapping;
} CDSArguments;

static void *gt_cds_arguments_new(void)
{
  CDSArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->seqfile = str_new();
  arguments->regionmapping = str_new();
  return arguments;
}

static void gt_cds_arguments_delete(void *tool_arguments)
{
  CDSArguments *arguments = tool_arguments;
  if (!arguments) return;
  str_delete(arguments->regionmapping);
  str_delete(arguments->seqfile);
  ma_free(arguments);
}

static OptionParser* gt_cds_option_parser_new(void *tool_arguments)
{
  CDSArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *option;
  assert(arguments);

  op = option_parser_new("[option ...] GFF3_file", "Add CDS features to exon "
                         "features given in GFF3_file.");

  /* -seqfile and -regionmapping */
  seqid2file_options(op, arguments->seqfile, arguments->regionmapping);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  option_parser_set_min_max_args(op, 1, 1);

  return op;
}

static int gt_cds_runner(UNUSED int argc, const char **argv, int parsed_args,
                         void *tool_arguments, Error *err)
{
  GenomeStream *gff3_in_stream, *cds_stream = NULL, *gff3_out_stream = NULL;
  GenomeNode *gn;
  CDSArguments *arguments = tool_arguments;
  RegionMapping *regionmapping;
  int had_err = 0;

  error_check(err);
  assert(arguments);

  /* create gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments->verbose);

  /* create region mapping */
  regionmapping = seqid2file_regionmapping_new(arguments->seqfile,
                                               arguments->regionmapping, err);
  if (!regionmapping)
    had_err = -1;

  /* create CDS stream */
  if (!had_err) {
    cds_stream = cds_stream_new(gff3_in_stream, regionmapping,
                                GT_CDS_SOURCE_TAG);
    if (!cds_stream)
      had_err = -1;
  }

  /* create gff3 output stream */
  /* XXX: replace NULL with proper outfile */
  if (!had_err)
    gff3_out_stream = gff3_out_stream_new(cds_stream, NULL);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, err)) &&
         gn) {
    genome_node_rec_delete(gn);
  }

  /* free */
  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(cds_stream);
  genome_stream_delete(gff3_in_stream);

  return had_err;
}

Tool *gt_cds(void)
{
  return tool_new(gt_cds_arguments_new,
                  gt_cds_arguments_delete,
                  gt_cds_option_parser_new,
                  NULL,
                  gt_cds_runner);
}
