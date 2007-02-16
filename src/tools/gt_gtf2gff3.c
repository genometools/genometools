/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(int *parsed_args, bool *be_tolerant, int argc,
                            char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  op = option_parser_new("[gtf_file]",
                         "Parse GTF2.2 file and show it as GFF3.");
  /* -tolerant */
  option = option_new_bool("tolerant", "be tolerant when parsing the GTF file",
                           be_tolerant, false);
  option_parser_add_option(op, option);
  /* parse */
  oprval = option_parser_parse_max_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, err);
  option_parser_free(op);
  return oprval;
}

int gt_gtf2gff3(int argc, char *argv[], Error *err)
{
  GenomeStream *gtf_in_stream,
                *gff3_out_stream;
  GenomeNode *gn;
  int parsed_args;
  bool be_tolerant;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &be_tolerant, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

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

  return 0;
}
