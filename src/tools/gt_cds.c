/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  unsigned int verbose;
} CDS_arguments;

static int parse_options(CDS_arguments *arguments, int argc,
                         char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *option;

  op = option_parser_new("[option ...] GFF3_file sequence_file",
                         "Add CDS features to exon features given in GFF3_file "
                         "(which refers to sequence_file).");

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* parse */
  parsed_args = option_parser_parse_min_max_args(op, argc, argv, versionfunc, 2,
                                                 2);
  option_parser_free(op);

  return parsed_args;
}

int gt_cds(int argc, char *argv[])
{
  Genome_stream *gff3_in_stream,
                *cds_stream,
                *gff3_out_stream;
  Genome_node *gn;
  CDS_arguments arguments;
  int parsed_args;

  /* option parsing */
  parsed_args = parse_options(&arguments, argc, argv);

  /* create gff3 input stream */
  assert(parsed_args < argc);
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose);

  /* create CDS stream */
  assert(parsed_args + 1 < argc);
  cds_stream = cds_stream_new(gff3_in_stream, argv[parsed_args + 1], "gt_cds");

  /* create gff3 output stream */
  gff3_out_stream = gff3_out_stream_new(cds_stream, stdout);

  /* pull the features through the stream and free them afterwards */
  while ((gn = genome_stream_next_tree(gff3_out_stream, NULL)))
    genome_node_rec_free(gn);

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(cds_stream);
  genome_stream_free(gff3_in_stream);

  return EXIT_SUCCESS;
}
