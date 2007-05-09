/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"
#include <libgtext/feature_index.h>

typedef struct {
  bool verbose;
  long offset;
  GenFile *outfp;
} Gff3_features_arguments;

static OPrval parse_options(int *parsed_args, Gff3_features_arguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  OutputFileInfo *ofi;
  Option  *option;
  OPrval oprval;
  env_error_check(env);

  /* init */
  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Work with and test GFF3 feature indexes.",
                         env);
  ofi = outputfileinfo_new(env);

  /* -offset */
  option = option_new_long("offset",
                           "transform all features by the given offset",
                           &arguments->offset, UNDEF_LONG, env);
  /* XXX: do not make this a ``normal option'' until the necessary envor checks
     have been added to range_offset() in range.c */
  option_is_development_option(option);
  option_parser_add_option(op, option, env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, ofi, env);

  /* parse options */
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);

  /* free */
  outputfileinfo_delete(ofi, env);
  option_parser_delete(op, env);

  return oprval;
}

int gt_gff3_features(int argc, const char **argv, Env *env)
{
  GenomeStream *sort_stream, *gff3_in_stream, *feature_stream = NULL;
  Gff3_features_arguments arguments;
  GenomeNode *gn;
  FeatureIndex *features = NULL;
  int parsed_args, has_err;
  unsigned long start, end;
  env_error_check(env);


  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments.verbose &&
                                               arguments.outfp, env);


  /* set offset (if necessary) */
  if (arguments.offset != UNDEF_LONG)
    gff3_in_stream_set_offset(gff3_in_stream, arguments.offset);
	
  /*  create sort stream */
  sort_stream = sort_stream_new(gff3_in_stream, env);
  
  /* create feature index and stream */
  features = feature_index_new(env);
  feature_stream = feature_stream_new(sort_stream, features, env);

  /* create gff3 output stream */
  /* gff3_out_stream = gff3_out_stream_new(feature_stream, arguments.outfp, env); */

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(feature_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn, env);
  }

  env_error_check(env);

  feature_index_print_contents(features, env); 

  /* free */
  feature_index_delete(features, env);
  genome_stream_delete(feature_stream, env);
  genome_stream_delete(sort_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  genfile_xclose(arguments.outfp, env);

  return has_err;
}
