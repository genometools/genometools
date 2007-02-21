/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool sort,
       mergefeat,
       verbose;
  long offset;
  FILE *outfp;
} Gff3_arguments;

static OPrval parse_options(int *parsed_args, Gff3_arguments *arguments,
                            int argc, char **argv, Env *env)
{
  OptionParser *op;
  Option *sort_option, *mergefeat_option, *option;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Parse, possibly transform, and output GFF3 files.");
  sort_option = option_new_bool("sort", "sort the GFF3 features (memory "
                                "consumption is O(file_size))",
                                &arguments->sort, false);
  option_parser_add_option(op, sort_option);
  mergefeat_option = option_new_bool("mergefeat", "merge adjacent features of "
                                     "the same type", &arguments->mergefeat,
                                     false);
  option_is_development_option(mergefeat_option);
  option_imply(mergefeat_option, sort_option);
  option_parser_add_option(op, mergefeat_option);
  option = option_new_long("offset",
                           "transform all features by the given offset",
                           &arguments->offset, UNDEFLONG);
  /* XXX: do not make this a ``normal option'' until the necessary envor checks
     have been added to range_offset() in range.c */
  option_is_development_option(option);
  option_parser_add_option(op, option);
  option = option_new_outputfile(&arguments->outfp);
  option_parser_add_option(op, option);
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op);
  return oprval;
}

int gt_gff3(int argc, char *argv[], Env *env)
{
  GenomeStream *gff3_in_stream,
                *sort_stream = NULL,
                *mergefeat_stream = NULL,
                *gff3_out_stream;
  Gff3_arguments arguments;
  GenomeNode *gn;
  int parsed_args, has_err;
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
                                               arguments.outfp != stdout);

  /* set offset (if necessary) */
  if (arguments.offset != UNDEFLONG)
    gff3_in_stream_set_offset(gff3_in_stream, arguments.offset);

  /* create sort stream (if necessary) */
  if (arguments.sort)
    sort_stream = sort_stream_new(gff3_in_stream);

  /* create merge feature stream (if necessary) */
  if (arguments.mergefeat) {
    assert(sort_stream);
    mergefeat_stream = mergefeat_stream_sorted_new(sort_stream);
  }

  /* create gff3 output stream */
  gff3_out_stream = gff3_out_stream_new(mergefeat_stream ? mergefeat_stream
                                        : (sort_stream ? sort_stream
                                                        : gff3_in_stream),
                                        arguments.outfp);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_out_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn);
  }

  /* free */
  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(sort_stream);
  genome_stream_delete(mergefeat_stream);
  genome_stream_delete(gff3_in_stream);
  if (arguments.outfp != stdout)
    xfclose(arguments.outfp);

  return has_err;
}
