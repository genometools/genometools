/*
  Copyright (c) Sascha Steinbiss, Malte Mader, Christin Schaerfer
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"
#include <libgtext/feature_index.h>

typedef struct {
  bool verbose;
  Str *seqid;
  unsigned long start, end;  
  GenFile *outfp;
} Gff3_view_arguments;

static OPrval parse_options(int *parsed_args, Gff3_view_arguments *arguments,
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

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);
  
  /* -seqid */
  arguments->seqid = str_new(env);
  option = option_new_string("seqid", "sequence region identifier", 
                            arguments->seqid,
                            "", env);
  option_parser_add_option(op, option, env);
  option_is_mandatory(option);
  option_hide_default(option);
  
  /* -start */
  option = option_new_ulong("start", "start position", 
                            &arguments->start,
                            1, env);
  option_parser_add_option(op, option, env);
  option_is_mandatory(option);
  option_hide_default(option);

  /* -end */
  option = option_new_ulong("end", "end position", 
                            &arguments->end,
                            1, env);
  option_parser_add_option(op, option, env);
  option_is_mandatory(option);
  option_hide_default(option);
  
  /* output file options */
  outputfile_register_options(op, &arguments->outfp, ofi, env);

  /* parse options */
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);

  /* free */
  outputfileinfo_delete(ofi, env);
  option_parser_delete(op, env);

  return oprval;
}

int gt_view(int argc, const char **argv, Env *env)
{
  GenomeStream *sort_stream, *gff3_in_stream, *feature_stream = NULL;
  Gff3_view_arguments arguments;
  GenomeNode *gn = NULL;
  FeatureIndex *features = NULL;
  int parsed_args, has_err;
  unsigned long i;
  Range qry_range;
  Array *results = NULL;
  env_error_check(env);


  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: 
      str_delete(arguments.seqid,env); 
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT: 
      str_delete(arguments.seqid,env); 
      return 0;
  }

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments.verbose &&
                                               arguments.outfp, env);

  /*  create sort stream */
  sort_stream = sort_stream_new(gff3_in_stream, env);
  
  /* create feature index and stream */
  features = feature_index_new(env);
  feature_stream = feature_stream_new(sort_stream, features, env);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(feature_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn, env);
  }

  /* sequence region id does not exist in gff file */
  if (!has_err && !feature_index_has_seqid(features, str_get(arguments.seqid), env)) {
    env_error_set(env, "sequence region '%s' does not exist in GFF file.", str_get(arguments.seqid));
    has_err = -1;
  }

  /* range end < range start */
  if (!has_err && (arguments.end < arguments.start)) {
    env_error_set(env, "end of query range precedes start of query range.");
    has_err = -1;
  }

  if(!has_err) {
    results = array_new(sizeof(GenomeNode*), env);
  
    qry_range.start = arguments.start;
    qry_range.end = arguments.end;
  
    feature_index_get_features_for_range(features, 
                                       results, 
                                       str_get(arguments.seqid), 
                                       qry_range, 
                                       env);

    printf("# of results: %lu\n", array_size(results)); 
  
    for(i=0;i<array_size(results);i++)
    {
      GenomeFeature *gf= **(GenomeFeature***) array_get(results, i);
      GenomeNode *gn = (GenomeNode*) gf;
      genome_node_traverse_children(gn,
                                  NULL,
                                  genome_node_print_feature_children,
                                  true,
                                  env); 
		  printf("------\n");
    }

    array_delete(results, env);  
	}
	
  /* free */
  str_delete(arguments.seqid,env);
  feature_index_delete(features, env);
  genome_stream_delete(feature_stream, env);
  genome_stream_delete(sort_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  genfile_xclose(arguments.outfp, env);

  return has_err;
}
