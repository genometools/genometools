/*
  Copyright (c) 2007-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/versionfunc.h"
#include "extended/chseqids_stream.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtdatahelp.h"
#include "extended/sort_stream_api.h"
#include "tools/gt_chseqids.h"

#define DEFAULT_JOINLENGTH 300

typedef struct {
  bool sort,
       verbose;
  GtFile *outfp;
} ChseqidsArguments;

static GtOPrval parse_options(int *parsed_args, ChseqidsArguments *arguments,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOutputFileInfo *ofi;
  GtOption *option;
  GtOPrval oprval;
  gt_error_check(err);

  /* init */
  op = gt_option_parser_new("[option ...] mapping_file [GFF3_file]",
                         "Change sequence ids by the mapping given in "
                         "mapping_file.");
  ofi = gt_output_file_info_new();

  /* -sort */
  option = gt_option_new_bool("sort",
                              "sort the GFF3 features after changing the "
                              "sequence ids\n(memory consumption is "
                              "proportional to the input file size)",
                              &arguments->sort, false);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, ofi);

  /* parse options */
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  gt_option_parser_set_min_max_args(op, 1, 2);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);

  /* free */
  gt_output_file_info_delete(ofi);
  gt_option_parser_delete(op);

  return oprval;
}

int gt_chseqids(int argc, const char **argv, GtError *err)
{
  GtNodeStream *gff3_in_stream, *chseqids_stream, *sort_stream = NULL,
               *gff3_out_stream = NULL;
  ChseqidsArguments arguments;
  GtStr *chseqids;
  int parsed_args, had_err = 0;

  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }

  /* create the streams */
  gff3_in_stream = gt_gff3_in_stream_new_sorted(argv[parsed_args + 1]);
  if (arguments.verbose && arguments.outfp)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);
  chseqids = gt_str_new_cstr(argv[parsed_args]);
  chseqids_stream = gt_chseqids_stream_new(gff3_in_stream, chseqids, err);
  if (!chseqids_stream)
    had_err = -1;
  gt_str_delete(chseqids);
  if (!had_err) {
    if (arguments.sort) {
      sort_stream = gt_sort_stream_new(chseqids_stream);
      gff3_out_stream = gt_gff3_out_stream_new(sort_stream, arguments.outfp);
    }
    else {
      gff3_out_stream = gt_gff3_out_stream_new(chseqids_stream,
                                               arguments.outfp);
    }
  }

  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(chseqids_stream);
  gt_node_stream_delete(sort_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_file_delete(arguments.outfp);

  return had_err;
}
