/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

#define GT_CDS_SOURCE_TAG "gt cds"

typedef struct {
  bool verbose;
} CDS_arguments;

static OPrval parse_options(int *parsed_args, CDS_arguments *arguments,
                            int argc, char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] GFF3_file sequence_file",
                         "Add CDS features to exon features given in GFF3_file "
                         "(refers to sequence_file).", env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* parse */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_cds(int argc, char *argv[], Env *env)
{
  GenomeStream *gff3_in_stream, *cds_stream, *gff3_out_stream = NULL;
  GenomeNode *gn;
  CDS_arguments arguments;
  int parsed_args, has_err = 0;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create gff3 input stream */
  assert(parsed_args < argc);
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose, env);

  /* create CDS stream */
  assert(parsed_args + 1 < argc);
  cds_stream = cds_stream_new(gff3_in_stream, argv[parsed_args + 1],
                              GT_CDS_SOURCE_TAG, env);
  if (!cds_stream)
    has_err = -1;

  /* create gff3 output stream */
  if (!has_err)
    gff3_out_stream = gff3_out_stream_new(cds_stream, stdout, env);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_out_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn, env);
  }

  /* free */
  genome_stream_delete(gff3_out_stream, env);
  genome_stream_delete(cds_stream, env);
  genome_stream_delete(gff3_in_stream, env);

  return has_err;
}
