/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool verbose;
} CDS_arguments;

static OPrval parse_options(int *parsed_args, CDS_arguments *arguments,
                            int argc, char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] GFF3_file sequence_file",
                         "Add CDS features to exon features given in GFF3_file "
                         "(which refers to sequence_file).");

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* parse */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, err);
  option_parser_free(op);
  return oprval;
}

int gt_cds(int argc, char *argv[], Error *err)
{
  GenomeStream *gff3_in_stream,
                *cds_stream,
                *gff3_out_stream;
  GenomeNode *gn;
  CDS_arguments arguments;
  int parsed_args, has_err;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

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
  while (!(has_err = genome_stream_next_tree(gff3_out_stream, &gn, NULL, err))
         && gn) {
    genome_node_rec_free(gn);
  }

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(cds_stream);
  genome_stream_free(gff3_in_stream);

  return has_err;
}
