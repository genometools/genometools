/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/option.h"
#include "core/outputfile.h"
#include "core/versionfunc.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/uniq_stream.h"
#include "tools/gt_uniq.h"

typedef struct {
  bool verbose;
  GtFile *outfp;
} UniqArguments;

static GtOPrval parse_options(int *parsed_args, UniqArguments *arguments,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOutputFileInfo *ofi;
  GtOption *option;
  GtOPrval oprval;
  gt_error_check(err);

  /* init */
  op = gt_option_parser_new("[option ...] [GFF3_file]", "Filter out repeated "
                         "features in a sorted GFF3_file.");
  ofi = gt_outputfileinfo_new();

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_outputfile_register_options(op, &arguments->outfp, ofi);

  /* parse options */
  gt_option_parser_set_max_args(op, 1);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);

  /* free */
  gt_outputfileinfo_delete(ofi);
  gt_option_parser_delete(op);

  return oprval;
}

int gt_uniq(int argc, const char **argv, GtError *err)
{
  GtNodeStream *gff3_in_stream,
               *uniq_stream = NULL,
               *gff3_out_stream = NULL;
  UniqArguments arguments;
  int parsed_args, had_err;
  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }

  /* create gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args]);
  if (arguments.verbose && arguments.outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create uniq stream */
  uniq_stream = gt_uniq_stream_new(gff3_in_stream);

  /* create gff3 output stream */
  gff3_out_stream = gt_gff3_out_stream_new(uniq_stream, arguments.outfp);

  /* pull the features through the stream and free them afterwards */
  had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(uniq_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_file_delete(arguments.outfp);

  return had_err;
}
