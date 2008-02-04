/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/stat_stream.h"
#include "tools/gt_stat.h"

typedef struct {
  bool verbose,
       gene_length_distribution,
       gene_score_distribution,
       exon_number_distribution,
       exon_length_distribution,
       intron_length_distribution;
} StatArguments;

static OPrval parse_options(int *parsed_args, StatArguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Show statistics about features contained in GFF3 "
                         "files.");

  /* -genelengthdistri */
  option = option_new_bool("genelengthdistri",
                           "show gene length distribution",
                           &arguments->gene_length_distribution, false);
  option_parser_add_option(op, option);

  /* -genescoresdistri */
  option = option_new_bool("genescoredistri", "show gene score distribution",
                           &arguments->gene_score_distribution, false);
  option_parser_add_option(op, option);

  /* -exonlengthdistri */
  option = option_new_bool("exonlengthdistri",
                           "show exon length distribution",
                           &arguments->exon_length_distribution, false);
  option_parser_add_option(op, option);

  /* -exonnumberdistri */
  option = option_new_bool("exonnumberdistri",
                           "show exon number distribution",
                           &arguments->exon_number_distribution, false);
  option_parser_add_option(op, option);

  /* -intronlengthdistri */
  option = option_new_bool("intronlengthdistri",
                           "show intron length distribution",
                           &arguments->intron_length_distribution, false);
  option_parser_add_option(op, option);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* parse */
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_stat(int argc, const char **argv, Error *err)
{
  GenomeStream *gff3_in_stream, *stat_stream;
  GenomeNode *gn;
  int parsed_args, had_err;
  StatArguments arguments;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments.verbose, false);

  /* create s status stream */
  stat_stream = stat_stream_new(gff3_in_stream,
                                arguments.gene_length_distribution,
                                arguments.gene_score_distribution,
                                arguments.exon_length_distribution,
                                arguments.exon_number_distribution,
                                arguments.intron_length_distribution);

  /* pull the features through the stream , compute the statistics, and free
     them afterwards */
  while (!(had_err = genome_stream_next_tree(stat_stream, &gn, err)) && gn) {
    genome_node_rec_delete(gn);
  }

  /* show statistics */
  if (!had_err)
    stat_stream_show_stats(stat_stream);

  /* free */
  genome_stream_delete(stat_stream);
  genome_stream_delete(gff3_in_stream);

  return had_err;
}
