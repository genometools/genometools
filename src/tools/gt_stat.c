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

#include "core/option.h"
#include "core/versionfunc.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"
#include "extended/stat_stream.h"
#include "tools/gt_stat.h"

typedef struct {
  bool verbose,
       gene_length_distribution,
       gene_score_distribution,
       exon_number_distribution,
       exon_length_distribution,
       intron_length_distribution;
} StatArguments;

static GtOPrval parse_options(int *parsed_args, StatArguments *arguments,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *option;
  GtOPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("[option ...] [GFF3_file ...]",
                         "Show statistics about features contained in GFF3 "
                         "files.");

  /* -genelengthdistri */
  option = gt_option_new_bool("genelengthdistri",
                           "show gene length distribution",
                           &arguments->gene_length_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -genescoresdistri */
  option = gt_option_new_bool("genescoredistri", "show gene score distribution",
                           &arguments->gene_score_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -exonlengthdistri */
  option = gt_option_new_bool("exonlengthdistri",
                           "show exon length distribution",
                           &arguments->exon_length_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -exonnumberdistri */
  option = gt_option_new_bool("exonnumberdistri",
                           "show exon number distribution",
                           &arguments->exon_number_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -intronlengthdistri */
  option = gt_option_new_bool("intronlengthdistri",
                           "show intron length distribution",
                           &arguments->intron_length_distribution, false);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* parse */
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_stat(int argc, const char **argv, GtError *err)
{
  GtNodeStream *gff3_in_stream, *stat_stream;
  int parsed_args, had_err;
  StatArguments arguments;
  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create a gff3 input stream */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                  argv + parsed_args);
  if (arguments.verbose)
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) gff3_in_stream);

  /* create s status stream */
  stat_stream = gt_stat_stream_new(gff3_in_stream,
                                   arguments.gene_length_distribution,
                                   arguments.gene_score_distribution,
                                   arguments.exon_length_distribution,
                                   arguments.exon_number_distribution,
                                   arguments.intron_length_distribution);

  /* pull the features through the stream , compute the statistics, and free
     them afterwards */
  had_err = gt_node_stream_pull(stat_stream, err);

  /* show statistics */
  if (!had_err)
    gt_stat_stream_show_stats(stat_stream);

  /* free */
  gt_node_stream_delete(stat_stream);
  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}
