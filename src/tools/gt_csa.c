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
  FILE *outfp;
} Csa_arguments;

static OPrval parse_options(int *parsed_args, Csa_arguments *arguments,
                            int argc, char *argv[], Error *err)
{
  OptionParser *op = option_parser_new("[option ...] [GFF3_file]",
                                       "Replace spliced alignments with "
                                       "computed consensus spliced "
                                       "alignments.");
  Option *option;
  OPrval oprval;
  error_check(err);

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
  oprval = option_parser_parse_max_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, err);
  option_parser_free(op);

  return oprval;
}

int gt_csa(int argc, char *argv[], Error *err)
{
  GenomeStream *gff3_in_stream,
                *csa_stream,
                *gff3_out_stream;
  GenomeNode *gn;
  Csa_arguments arguments;
  Log *log = NULL;
  int parsed_args, has_err;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create log object if necessary */
  if (arguments.debug)
    log = log_new();

  /* create the streams */
  gff3_in_stream  = gff3_in_stream_new_sorted(argv[parsed_args],
                                              arguments.verbose &&
                                              arguments.outfp != stdout);
  csa_stream      = csa_stream_new(gff3_in_stream, arguments.join_length);
  gff3_out_stream = gff3_out_stream_new(csa_stream, arguments.outfp);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_out_stream, &gn, log, err))
         && gn) {
    genome_node_rec_free(gn);
  }

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(csa_stream);
  genome_stream_free(gff3_in_stream);
  log_free(log);
  if (arguments.outfp != stdout)
    xfclose(arguments.outfp);

  return has_err;
}
