/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <gt.h>

typedef struct {
  Str *sequence_region_id;
  unsigned long from,
                to;
  int width;
  bool pipe,
       verbose;
} Png_info;

static OPrval parse_options(int *parsed_args, Png_info *info, int argc,
                            const char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  bool force;
  env_error_check(env);
  op = option_parser_new("[option ...] png_file [gff3_file]",
                         "Draw gff3 features as png.", env);

  /* -sequence-region */
  option = option_new_string("sequence-region", "set drawed sequence-region "
                             "(by default the first one is taken)",
                             info->sequence_region_id, NULL, env);
  option_hide_default(option);
  option_parser_add_option(op, option, env);

  /* -from */
  option = option_new_ulong_min("from", "set position from where on features "
                                "are shown", &info->from, 1, 1, env);
  option_hide_default(option);
  option_parser_add_option(op, option, env);

  /* -to */
  option = option_new_ulong("to", "set maximal position of features to be "
                            "shown", &info->to, ~0UL, env);
  option_hide_default(option);
  option_parser_add_option(op, option, env);

  /* -width */
  option = option_new_int_min("width", "set the width of the png file",
                              &info->width, 1024, 1, env);
  option_parser_add_option(op, option, env);

  /* -pipe */
  option = option_new_bool("pipe", "use pipe mode (i.e., show all gff3 "
                           "features on stdout)", &info->pipe, false, env);
  option_parser_add_option(op, option, env);

  /* -force */
  option = option_new_bool("force", "force writing to output file", &force,
                           false, env);
  option_parser_add_option(op, option, env);

  /* -v */
  option = option_new_verbose(&info->verbose, env);
  option_parser_add_option(op, option, env);

  /* parsing */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 2, env);

  /* checks */
  if (oprval == OPTIONPARSER_OK && info->from > info->to) {
    env_error_set(env, "\"-from\" must be <= \"-to\"");
    oprval = OPTIONPARSER_ERROR;
  }

  if (oprval == OPTIONPARSER_OK && !force && file_exists(argv[*parsed_args])) {
    env_error_set(env, "file \"%s\" exists already. use option -force to "
                  "overwrite", argv[*parsed_args]);
    oprval = OPTIONPARSER_ERROR;
  }

  option_parser_delete(op, env);
  return oprval;
}

int gt_png(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream, *png_stream, *gff3_out_stream = NULL;
  GenomeNode *gn;
  int parsed_args, has_err;
  Png_info info;

  /* option parsing */
  info.sequence_region_id = str_new(env);
  switch (parse_options(&parsed_args, &info, argc, argv, env)) {
    case OPTIONPARSER_OK:
      break;
    case OPTIONPARSER_ERROR:
      str_delete(info.sequence_region_id, env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(info.sequence_region_id, env);
      return 0;
  }

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args + 1],
                                             info.verbose && !info.pipe, env);

  /* create a cairo stream */
  png_stream = png_stream_new(gff3_in_stream, info.sequence_region_id,
                                  info.from, info.to, argv[parsed_args],
                                  info.width, env);

  /* create gff3 output stream if -pipe was used */
  if (info.pipe)
    gff3_out_stream = gff3_out_stream_new(png_stream, NULL, env);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(info.pipe ? gff3_out_stream
                                                       : png_stream,
                                             &gn, env)) && gn) {
    genome_node_rec_delete(gn, env);
  }

  /* draw */
  if (!has_err)
    png_stream_draw((PNGStream*) png_stream, info.verbose, env);

  /* free */
  genome_stream_delete(gff3_out_stream, env);
  genome_stream_delete(png_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  str_delete(info.sequence_region_id, env);

  return has_err;
}
