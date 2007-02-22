/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool verbose,
       gene_length_distribution,
       gene_score_distribution,
       exon_length_distribution,
       intron_length_distribution;
} Stat_arguments;

typedef struct {
  unsigned long number_of_trees;
  GenomeVisitor *stat_visitor;
} StatInfo;

static OPrval parse_options(int *parsed_args, Stat_arguments *arguments,
                            int argc, char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Show statistics about features contained in GFF3 "
                         "files.", env);

  /* -genelengthdistri */
  option = option_new_bool("genelengthdistri",
                           "show gene length distribution",
                           &arguments->gene_length_distribution, false, env);
  option_parser_add_option(op, option, env);

  /* -genescoresdistri */
  option = option_new_bool("genescoredistri", "show gene score distribution",
                           &arguments->gene_score_distribution, false, env);
  option_parser_add_option(op, option, env);

  /* -exonlengthdistri */
  option = option_new_bool("exonlengthdistri",
                           "show exon length distribution",
                           &arguments->exon_length_distribution, false, env);
  option_parser_add_option(op, option, env);

  /* -intronlengthdistri */
  option = option_new_bool("intronlengthdistri",
                           "show intron length distribution",
                           &arguments->intron_length_distribution, false, env);
  option_parser_add_option(op, option, env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* parse */
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

static int compute_statistics(GenomeNode *gn, void *data, Env *env)
{
  StatInfo *info = (StatInfo*) data;
  env_error_check(env);
  assert(info && info->stat_visitor);
  return genome_node_accept(gn, info->stat_visitor, env);
}

int gt_stat(int argc, char *argv[], Env *env)
{
  GenomeStream *gff3_in_stream;
  GenomeNode *gn;
  int has_err, parsed_args;
  Stat_arguments arguments;
  StatInfo info;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* init */
  info.number_of_trees = 0;
  info.stat_visitor = stat_visitor_new(arguments.gene_length_distribution,
                                       arguments.gene_score_distribution,
                                       arguments.exon_length_distribution,
                                       arguments.intron_length_distribution);

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments.verbose, env);

  /* pull the features through the stream , compute the statistics, and free
     them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_in_stream, &gn, env)) && gn) {
    info.number_of_trees++;
    has_err = genome_node_traverse_children(gn, &info, compute_statistics, true,
                                            env);
    genome_node_rec_delete(gn, env);
    if (has_err)
      break;
  }

  /* show statistics */
  if (!has_err) {
    printf("parsed feature trees: %lu\n", info.number_of_trees);
    stat_visitor_show_stats(info.stat_visitor);
  }

  /* free */
  genome_visitor_delete(info.stat_visitor, env);
  genome_stream_delete(gff3_in_stream, env);

  return has_err;
}
