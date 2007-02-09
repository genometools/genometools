/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  unsigned int join,
               translate,
               verbose;
  Str *type;
} Extractfeat_arguments;

static int parse_options(Extractfeat_arguments *arguments, int argc,
                         char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *option;

  op = option_parser_new("[option ...] GFF3_file sequence_file",
                         "Extract features given in GFF3_file from "
                         "sequence_file.");

  /* -type */
  option = option_new_string("type", "set type of features to extract",
                             arguments->type, NULL);
  option_is_mandatory(option);
  option_parser_add_option(op, option);

  /* -join */
  option = option_new_boolean("join", "join feature sequences in the same "
                              "subgraph into a single one", &arguments->join,
                              0);
  option_parser_add_option(op, option);

  /* -translate */
  option = option_new_boolean("translate", "translate the features (of a DNA "
                              "sequence) into protein", &arguments->translate,
                              0);
  option_parser_add_option(op, option);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* parse */
  parsed_args = option_parser_parse_min_max_args(op, argc, argv, versionfunc,
                                                 2, 2);
  option_parser_free(op);

  return parsed_args;
}

int gt_extractfeat(int argc, char *argv[])
{
  Genome_stream *gff3_in_stream,
                *extractfeat_stream;
  Genome_node *gn;
  Genome_feature_type type;
  Extractfeat_arguments arguments;
  int parsed_args;

  /* option parsing */
  arguments.type = str_new();
  parsed_args = parse_options(&arguments, argc, argv);

  /* determine type and make sure it is a valid one */
  if (genome_feature_type_get(&type, str_get(arguments.type)) == -1)
    error("\"%s\" is not a valid feature type", str_get(arguments.type));

  /* create gff3 input stream */
  assert(parsed_args < argc);
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose);

  /* create extract feature stream */
  assert(parsed_args + 1 < argc);
  extractfeat_stream = extractfeat_stream_new(gff3_in_stream,
                                              argv[parsed_args + 1],
                                              type, arguments.join,
                                              arguments.translate);

  /* pull the features through the stream and free them afterwards */
  while ((gn = genome_stream_next_tree(extractfeat_stream, NULL)))
    genome_node_rec_free(gn);

  /* free */
  genome_stream_free(extractfeat_stream);
  genome_stream_free(gff3_in_stream);
  str_free(arguments.type);

  return EXIT_SUCCESS;
}
