/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  unsigned int verbose;
  unsigned long max_gene_length;
  double min_gene_score;
  FILE *outfp;
} Filter_arguments;

static int parse_options(Filter_arguments *arguments, int argc, char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *option;

  op = option_parser_new("[option ...] [GFF3_file ...]", "Filter GFF3 files.");

  /* -maxgenelength */
  option = option_new_ulong_min("maxgenelength", "the maximum length a gene "
                                "can have to pass the filter",
                                &arguments->max_gene_length, UNDEFULONG, 1);
  option_parser_add_option(op, option);

  /* -mingenescore */
  option = option_new_double("mingenescore", "the minimum score a gene must "
                             "have to pass the filter",
                             &arguments->min_gene_score, UNDEFDOUBLE);
  option_parser_add_option(op, option);

  /* -o */
  option = option_new_outputfile(&arguments->outfp);
  option_parser_add_option(op, option);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* parse */
  parsed_args = option_parser_parse(op, argc, argv, versionfunc);
  option_parser_free(op);

  return parsed_args;
}

int gt_filter(int argc, char *argv[])
{
  Genome_stream *gff3_in_stream,
                *filter_stream,
                *gff3_out_stream;
  Genome_node *gn;
  Filter_arguments arguments;
  int parsed_args;

  /* option parsing */
  parsed_args = parse_options(&arguments, argc, argv);

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments.verbose &&
                                               arguments.outfp != stdout);

  /* create a filter stream */
  filter_stream = filter_stream_new(gff3_in_stream, arguments.max_gene_length,
                                    arguments.min_gene_score);

  /* create a gff3 output stream */
  gff3_out_stream = gff3_out_stream_new(filter_stream, arguments.outfp);

  /* pull the features through the stream and free them afterwards */
  while ((gn = genome_stream_next_tree(gff3_out_stream, NULL)))
    genome_node_rec_free(gn);

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(filter_stream);
  genome_stream_free(gff3_in_stream);
  if (arguments.outfp != stdout)
    xfclose(arguments.outfp);

  return EXIT_SUCCESS;
}
