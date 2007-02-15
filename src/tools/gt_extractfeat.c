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

static int parse_options(ExtractFeatArguments *arguments, int argc, char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *option, *seqfile_option, *regionmapping_option;

  op = option_parser_new("[option ...] GFF3_file",
                         "Extract features given in GFF3_file from "
                         "sequence file.");

  /* -type */
  option = option_new_string("type", "set type of features to extract",
                             arguments->type, NULL);
  option_is_mandatory(option);
  option_parser_add_option(op, option);

  /* -join */
  option = option_new_bool("join", "join feature sequences in the same "
                           "subgraph into a single one", &arguments->join,
                           false);
  option_parser_add_option(op, option);

  /* -translate */
  option = option_new_bool("translate", "translate the features (of a DNA "
                           "sequence) into protein", &arguments->translate,
                           false);
  option_parser_add_option(op, option);

  /* -seqfile */
  seqfile_option = option_new_string("seqfile", "set the sequence file from "
                                     "which to extract the features",
                                     arguments->seqfile, NULL);
  option_parser_add_option(op, seqfile_option);

  /* -regionmapping */
  regionmapping_option = option_new_string("regionmapping", "set file "
                                           "containing sequence-region to "
                                           "sequence file mapping",
                                           arguments->regionmapping, NULL);
  option_parser_add_option(op, regionmapping_option);

  /* either option -seqfile or -regionmapping is mandatory */
  option_is_mandatory_either(seqfile_option, regionmapping_option);

  /* the options -seqfile and -regionmapping exclude each other */
  option_exclude(seqfile_option, regionmapping_option);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  parsed_args = option_parser_parse_min_max_args(op, argc, argv, versionfunc,
                                                 1, 1);
  option_parser_free(op);

  return parsed_args;
}

int gt_extractfeat(int argc, char *argv[])
{
  GenomeStream *gff3_in_stream,
                *extractfeat_stream;
  GenomeNode *gn;
  GenomeFeatureType type;
  ExtractFeatArguments arguments;
  RegionMapping *regionmapping;
  int parsed_args;
  Error *err = error_new();
  int has_err = 0;

  /* option parsing */
  arguments.type = str_new();
  arguments.seqfile = str_new();
  arguments.regionmapping = str_new();
  parsed_args = parse_options(&arguments, argc, argv);

  /* determine type and make sure it is a valid one */
  if (genome_feature_type_get(&type, str_get(arguments.type)) == -1)
    error("\"%s\" is not a valid feature type", str_get(arguments.type));

  /* create gff3 input stream */
  assert(parsed_args < argc);
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose);

  /* create extract feature stream */
  extractfeat_stream = extractfeat_stream_new(gff3_in_stream,
                                              type, arguments.join,
                                              arguments.translate);

  /* set sequence source */
  assert(str_get(arguments.seqfile) || str_get(arguments.regionmapping));
  assert(!(str_get(arguments.seqfile) && str_get(arguments.regionmapping)));
  if (str_get(arguments.seqfile))
    extractfeat_stream_use_sequence_file(extractfeat_stream, arguments.seqfile);
  else {
    regionmapping = regionmapping_new(arguments.regionmapping, err);
    if (!regionmapping)
      has_err = -1;
    else
      extractfeat_stream_use_region_mapping(extractfeat_stream, regionmapping);
  }

  /* pull the features through the stream and free them afterwards */
  if (!has_err) {
    while (!genome_stream_next_tree(extractfeat_stream, &gn, NULL, err) && gn)
      genome_node_rec_free(gn);
  }

  /* free */
  genome_stream_free(extractfeat_stream);
  genome_stream_free(gff3_in_stream);
  str_free(arguments.regionmapping);
  str_free(arguments.seqfile);
  str_free(arguments.type);
  error_abort(err);
  error_free(err);

  return EXIT_SUCCESS;
}
