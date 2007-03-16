/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

#define DEFAULT_JOINLENGTH 300

typedef struct {
  bool verbose;
  unsigned long join_length;
  FILE *outfp;
} Csa_arguments;

static OPrval parse_options(int *parsed_args, Csa_arguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op = option_parser_new("[option ...] [GFF3_file]",
                                       "Replace spliced alignments with "
                                       "computed consensus spliced "
                                       "alignments.", env);
  Option *option;
  OPrval oprval;
  env_error_check(env);

  /* -join-length */
  option = option_new_ulong("join-length", "set join length for the spliced "
                            "alignment clustering", &arguments->join_length,
                            DEFAULT_JOINLENGTH, env);
  option_parser_add_option(op, option, env);

  /* -o */
  option = option_new_outputfile(&arguments->outfp, env);
  option_parser_add_option(op, option, env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* parse */
  oprval = option_parser_parse_max_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, env);
  option_parser_delete(op, env);

  return oprval;
}

int gt_csa(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream,
                *csa_stream,
                *gff3_out_stream;
  GenomeNode *gn;
  Csa_arguments arguments;
  int parsed_args, has_err;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create the streams */
  gff3_in_stream  = gff3_in_stream_new_sorted(argv[parsed_args],
                                              arguments.verbose &&
                                              arguments.outfp != stdout, env);
  csa_stream      = csa_stream_new(gff3_in_stream, arguments.join_length, env);
  gff3_out_stream = gff3_out_stream_new(csa_stream, arguments.outfp, env);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_out_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn, env);
  }

  /* free */
  genome_stream_delete(gff3_out_stream, env);
  genome_stream_delete(csa_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  if (arguments.outfp != stdout)
    xfclose(arguments.outfp);

  return has_err;
}
