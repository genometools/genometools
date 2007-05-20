/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  unsigned long max_feature_dist;
  bool verbose;
} RegionCovArguments;

static OPrval parse_options(int *parsed_args, RegionCovArguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] GFF3_file",
                         "Show which parts of the given sequence regions are "
                         "covered by features.", env);
  /* -maxfeaturedist */
  o = option_new_ulong("maxfeaturedist", "set the maximum distance two "
                       "features can have while still being in the same "
                       "``cluster''", &arguments->max_feature_dist, 0, env);
  option_parser_add_option(op, o, env);
  /* -v */
  o = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, o, env);
  /* parse */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_regioncov(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream;
  GenomeNode *gn;
  RegionCovArguments arguments;
  int parsed_args, has_err = 0;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      return 0;
  }

  /* create gff3 input stream */
  assert(parsed_args < argc);
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose, env);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_in_stream, &gn, env)) && gn) {
      genome_node_rec_delete(gn, env);
  }

  /* free */
  genome_stream_delete(gff3_in_stream, env);

  return has_err;
}
