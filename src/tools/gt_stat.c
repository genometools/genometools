/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool verbose,
       gene_length_distribution,
       gene_score_distribution;
} Stat_arguments;

typedef struct {
  unsigned long number_of_trees;
  GenomeVisitor *stat_visitor;
} StatInfo;

static OPrval parse_options(int *parsed_args, Stat_arguments *arguments,
                            int argc, char **argv, Error *err)
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

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* parse */
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_free(op);
  return oprval;
}

static int compute_statistics(GenomeNode *gn, void *data, Error *err)
{
  StatInfo *info = (StatInfo*) data;
  error_check(err);
  assert(info && info->stat_visitor);
  return genome_node_accept(gn, info->stat_visitor, NULL, err);
}

int gt_stat(int argc, char *argv[], Error *err)
{
  GenomeStream *gff3_in_stream;
  GenomeNode *gn;
  int has_err, parsed_args;
  Stat_arguments arguments;
  StatInfo info;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* init */
  info.number_of_trees = 0;
  info.stat_visitor = stat_visitor_new(arguments.gene_length_distribution,
                                       arguments.gene_score_distribution);

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments.verbose);

  /* pull the features through the stream , compute the statistics, and free
     them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_in_stream, &gn, NULL, err)) &&
         gn) {
    info.number_of_trees++;
    has_err = genome_node_traverse_children(gn, &info, compute_statistics, true,
                                            err);
    genome_node_rec_free(gn);
    if (has_err)
      break;
  }

  /* show statistics */
  if (!has_err) {
    printf("parsed feature trees: %lu\n", info.number_of_trees);
    stat_visitor_show_stats(info.stat_visitor);
  }

  /* free */
  genome_visitor_free(info.stat_visitor);
  genome_stream_free(gff3_in_stream);

  return has_err;
}
