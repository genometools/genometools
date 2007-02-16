/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

#define DEFAULT_JOINLENGTH 300

typedef struct {
  bool verbose,
       debug;
  unsigned long join_length;
  Log *log;
  FILE *outfp;
} Csa_arguments;

static int parse_options(Csa_arguments *arguments, int argc, char *argv[])
{
  int parsed_args;
  OptionParser *op = option_parser_new("[option ...] [GFF3_file]",
                                       "Replace spliced alignments with "
                                       "computed consensus spliced "
                                       "alignments.");
  Option *option;

  /* -join-length */
  option = option_new_ulong("join-length", "set join length for the spliced "
                            "alignment clustering", &arguments->join_length,
                            DEFAULT_JOINLENGTH);
  option_parser_add_option(op, option);

  /* -o */
  option = option_new_outputfile(&arguments->outfp);
  option_parser_add_option(op, option);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* -debug */
  option = option_new_debug(&arguments->debug);
  option_parser_add_option(op, option);

  /* parse */
  option_parser_parse_max_args(op, &parsed_args, argc, argv, versionfunc, 1);
  option_parser_free(op);

  if (arguments->debug)
    arguments->log = log_new();
  else
    arguments->log = NULL;

  return parsed_args;
}

int gt_csa(int argc, char *argv[])
{
  GenomeStream *gff3_in_stream,
                *csa_stream,
                *gff3_out_stream;
  GenomeNode *gn;
  Csa_arguments arguments;
  int parsed_args;
  Error *err = error_new();

  /* option parsing */
  parsed_args = parse_options(&arguments, argc, argv);

  /* create the streams */
  gff3_in_stream  = gff3_in_stream_new_sorted(argv[parsed_args],
                                              arguments.verbose &&
                                              arguments.outfp != stdout);
  csa_stream      = csa_stream_new(gff3_in_stream, arguments.join_length);
  gff3_out_stream = gff3_out_stream_new(csa_stream, arguments.outfp);

  /* pull the features through the stream and free them afterwards */
  while (!genome_stream_next_tree(gff3_out_stream, &gn, arguments.log, err) &&
         gn) {
    genome_node_rec_free(gn);
  }

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(csa_stream);
  genome_stream_free(gff3_in_stream);
  log_free(arguments.log);
  if (arguments.outfp != stdout)
    xfclose(arguments.outfp);

  error_abort(err);
  error_free(err);

  return EXIT_SUCCESS;
}
