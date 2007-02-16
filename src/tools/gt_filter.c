/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool verbose;
  unsigned long max_gene_length;
  double min_gene_score;
  FILE *outfp;
} FilterArgumentss;

static OPrval parse_options(int *parsed_args, FilterArgumentss *arguments,
                            int argc, char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  error_check(err);

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
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_free(op);
  return oprval;
}

int gt_filter(int argc, char *argv[], Error *err)
{
  GenomeStream *gff3_in_stream,
                *filter_stream,
                *gff3_out_stream;
  GenomeNode *gn;
  FilterArgumentss arguments;
  int parsed_args;

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

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
  while (!genome_stream_next_tree(gff3_out_stream, &gn, NULL, err) && gn)
    genome_node_rec_free(gn);

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(filter_stream);
  genome_stream_free(gff3_in_stream);
  if (arguments.outfp != stdout)
    xfclose(arguments.outfp);

  return 0;
}
