/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool join,
       translate,
       verbose;
  Str *type,
      *seqfile,
      *regionmapping;
} ExtractFeatArguments;

static OPrval parse_options(int *parsed_args, ExtractFeatArguments *arguments,
                            int argc, char **argv, Env *env)
{
  OptionParser *op;
  Option *option, *seqfile_option, *regionmapping_option;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] GFF3_file",
                         "Extract features given in GFF3_file from "
                         "sequence file.", env);

  /* -type */
  option = option_new_string("type", "set type of features to extract",
                             arguments->type, NULL, env);
  option_is_mandatory(option);
  option_parser_add_option(op, option, env);

  /* -join */
  option = option_new_bool("join", "join feature sequences in the same "
                           "subgraph into a single one", &arguments->join,
                           false, env);
  option_parser_add_option(op, option, env);

  /* -translate */
  option = option_new_bool("translate", "translate the features (of a DNA "
                           "sequence) into protein", &arguments->translate,
                           false, env);
  option_parser_add_option(op, option, env);

  /* -seqfile */
  seqfile_option = option_new_string("seqfile", "set the sequence file from "
                                     "which to extract the features",
                                     arguments->seqfile, NULL, env);
  option_parser_add_option(op, seqfile_option, env);

  /* -regionmapping */
  regionmapping_option = option_new_string("regionmapping", "set file "
                                           "containing sequence-region to "
                                           "sequence file mapping",
                                           arguments->regionmapping, NULL, env);
  option_parser_add_option(op, regionmapping_option, env);

  /* either option -seqfile or -regionmapping is mandatory */
  option_is_mandatory_either(seqfile_option, regionmapping_option);

  /* the options -seqfile and -regionmapping exclude each other */
  option_exclude(seqfile_option, regionmapping_option, env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_extractfeat(int argc, char *argv[], Env *env)
{
  GenomeStream *gff3_in_stream = NULL, *extractfeat_stream = NULL;
  GenomeNode *gn;
  GenomeFeatureType type;
  ExtractFeatArguments arguments;
  RegionMapping *regionmapping;
  int parsed_args;
  int has_err = 0;
  env_error_check(env);

  /* option parsing */
  arguments.type = str_new(env);
  arguments.seqfile = str_new(env);
  arguments.regionmapping = str_new(env);
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(arguments.regionmapping, env);
      str_delete(arguments.seqfile, env);
      str_delete(arguments.type, env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(arguments.regionmapping, env);
      str_delete(arguments.seqfile, env);
      str_delete(arguments.type, env);
      return 0;
  }

  /* determine type and make sure it is a valid one */
  if (genome_feature_type_get(&type, str_get(arguments.type)) == -1) {
    env_error_set(env, "\"%s\" is not a valid feature type",
                  str_get(arguments.type));
    has_err = -1;
  }

  if (!has_err) {
    /* create gff3 input stream */
    assert(parsed_args < argc);
    gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                               arguments.verbose, env);

    /* set sequence source */
    assert(str_length(arguments.seqfile) ||
           str_length(arguments.regionmapping));
    assert(!(str_length(arguments.seqfile) &&
             str_length(arguments.regionmapping)));
    if (str_length(arguments.seqfile))
      regionmapping = regionmapping_new_seqfile(arguments.seqfile, env);
    else
      regionmapping = regionmapping_new_mapping(arguments.regionmapping, env);
    if (!regionmapping)
      has_err = -1;
  }

  if (!has_err) {
    /* create extract feature stream */
    extractfeat_stream = extractfeat_stream_new(gff3_in_stream, regionmapping,
                                                type, arguments.join,
                                                arguments.translate, env);

    /* pull the features through the stream and free them afterwards */
    while (!(has_err = genome_stream_next_tree(extractfeat_stream, &gn, env)) &&
           gn) {
      genome_node_rec_delete(gn, env);
    }
  }

  /* free */
  genome_stream_delete(extractfeat_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  str_delete(arguments.regionmapping, env);
  str_delete(arguments.seqfile, env);
  str_delete(arguments.type, env);

  return has_err;
}
