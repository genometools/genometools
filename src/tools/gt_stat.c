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
  Genome_visitor *stat_visitor;
} Stat_info;

static int parse_options(Stat_arguments *arguments, int argc, char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *option;

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
  parsed_args = option_parser_parse(op, argc, argv, versionfunc);
  option_parser_free(op);

  return parsed_args;
}

static void compute_statistics(GenomeNode *gn, void *data)
{
  Stat_info *info = (Stat_info*) data;
  assert(info && info->stat_visitor);
  genome_node_accept(gn, info->stat_visitor, NULL);
}

int gt_stat(int argc, char *argv[])
{
  Genome_stream *gff3_in_stream;
  GenomeNode *gn;
  int parsed_args;
  Stat_arguments arguments;
  Stat_info info;

  /* option parsing */
  parsed_args = parse_options(&arguments, argc, argv);

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
  while ((gn = genome_stream_next_tree(gff3_in_stream, NULL))) {
    info.number_of_trees++;
    genome_node_traverse_children(gn, &info, compute_statistics, 1);
    genome_node_rec_free(gn);
  }

  /* show statistics */
  printf("parsed feature trees: %lu\n", info.number_of_trees);
  stat_visitor_show_stats(info.stat_visitor);

  /* free */
  genome_visitor_free(info.stat_visitor);
  genome_stream_free(gff3_in_stream);

  return EXIT_SUCCESS;
}
