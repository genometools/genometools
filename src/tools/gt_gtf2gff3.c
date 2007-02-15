/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static int parse_options(bool *be_tolerant, int argc, char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *option;

  op = option_parser_new("[gtf_file]",
                         "Parse GTF2.2 file and show it as GFF3.");

  /* -tolerant */
  option = option_new_bool("tolerant", "be tolerant when parsing the GTF file",
                           be_tolerant, false);
  option_parser_add_option(op, option);

  /* parse */
  parsed_args = option_parser_parse_max_args(op, argc, argv, versionfunc, 1);
  option_parser_free(op);

  return parsed_args;
}

int gt_gtf2gff3(int argc, char *argv[])
{
  GenomeStream *gtf_in_stream,
                *gff3_out_stream;
  GenomeNode *gn;
  int parsed_args;
  bool be_tolerant;
  Error *err = error_new();

  /* option parsing */
  parsed_args = parse_options(&be_tolerant, argc, argv);

  /* create a gtf input stream */
  gtf_in_stream = gtf_in_stream_new(argv[parsed_args], be_tolerant, err);
  if (!gtf_in_stream)
    error_abort(err);

  /* create a gff3 output stream */
  gff3_out_stream = gff3_out_stream_new(gtf_in_stream, stdout);

  /* pull the features through the stream and free them afterwards */
  while (!genome_stream_next_tree(gff3_out_stream, &gn, NULL, err) && gn)
    genome_node_rec_free(gn);

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(gtf_in_stream);
  error_abort(err);
  error_free(err);

  return EXIT_SUCCESS;
}
