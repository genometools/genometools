/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Malte Mader <mmader@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool verbose;
  Str *seqid, *outfile;
  unsigned long start, end;
  unsigned int width;
} Gff3_view_arguments;

static OPrval parse_options(int *parsed_args, Gff3_view_arguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option  *option;
  OPrval oprval;
  env_error_check(env);

  /* init */
  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Create visual representations of GFF annotations.",
                         env);
                         
  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* -o */
  arguments->outfile = str_new(env);
  option = option_new_string("o", "output file",
                            arguments->outfile,
                            "out.png", env);
  option_parser_add_option(op, option, env);
  option_is_mandatory(option);

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

	/* -width */
  option = option_new_uint("width", "image width",
                            &arguments->width,
                            600, env);
  option_parser_add_option(op, option, env);

  /* parse options */
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);

  /* free */
  option_parser_delete(op, env);

  return oprval;
}

int gt_view(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream = NULL,
               *feature_stream = NULL;
  Gff3_view_arguments arguments;
  GenomeNode *gn = NULL;
  FeatureIndex *features = NULL;
  int parsed_args, has_err=0;
  unsigned long i;
  Range qry_range;
  Array *results = NULL;
  Config *cfg;
  Str *config_file;
  
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(arguments.outfile,env);
      str_delete(arguments.seqid,env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(arguments.outfile,env);
      str_delete(arguments.seqid,env);
      return 0;
  }

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_sorted(*(argv + parsed_args),
                                             arguments.verbose &&
                                             NULL, env);

  /* create feature index and stream */
  features = feature_index_new(env);
  feature_stream = feature_stream_new(gff3_in_stream, features, env);

  /* check for correct order: range end < range start */
  if (!has_err && (arguments.end < arguments.start))
  {
    env_error_set(env, "end of query range (%lu) precedes "
		                   "start of query range (%lu)",
		                   arguments.end, arguments.start);
    has_err = -1;
  }

  if(!has_err)
  {
    /* pull the features through the stream and free them afterwards */
    while (!(has_err = genome_stream_next_tree(feature_stream, &gn, env)) && gn)
    {
      genome_node_rec_delete(gn, env);
    }
  }

  /* sequence region id does not exist in gff file */
  if (!has_err && !feature_index_has_seqid(features,
                                           str_get(arguments.seqid),
                                           env))
  {
    env_error_set(env, "sequence region '%s' does not exist in GFF input file",
		              str_get(arguments.seqid));
    has_err = -1;
  }

  if (!has_err) {
    results = array_new(sizeof (GenomeNode*), env);

    qry_range.start = arguments.start;
    qry_range.end = arguments.end;

    feature_index_get_features_for_range(features,
                                       results,
                                       str_get(arguments.seqid),
                                       qry_range,
                                       env);

    /* Silence is golden. */
    if (arguments.verbose)
    {
      printf("# of results: %lu\n", array_size(results));
      for (i=0;i<array_size(results);i++)
      {
        GenomeFeature *gf= *(GenomeFeature**) array_get(results, i);
        GenomeNode *gn = (GenomeNode*) gf;
        genome_node_traverse_children(gn,
                                      NULL,
                                      genome_node_print_feature_children,
                                      true,
                                      env);
	  	  printf("------\n");
      }
    }

    config_file = str_new_cstr("config.lua", env);
    
    cfg = config_new(env, arguments.verbose);
    config_load_file(cfg, config_file, env);

    Diagram* d = diagram_new(results, qry_range, cfg, env);
    Render* r = render_new(d, cfg, env);

    render_to_png(r, str_get(arguments.outfile), arguments.width, env);

    render_delete(r, env);
    config_delete(cfg, env);
    str_delete(config_file, env);
    diagram_delete(d, env);
    array_delete(results, env);
  }

  /* free */
  str_delete(arguments.seqid,env);
	str_delete(arguments.outfile,env);
  feature_index_delete(features, env);
  genome_stream_delete(feature_stream, env);
  genome_stream_delete(gff3_in_stream, env);

  return has_err;
}
