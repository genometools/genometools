/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Malte Mader <mmader@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool pipe,
       verbose;
  Str *seqid;
  unsigned long start,
                end;
  unsigned int width;
} Gff3_view_arguments;

static OPrval parse_options(int *parsed_args, Gff3_view_arguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option  *option, *option2;
  OPrval oprval;
  bool force;
  env_error_check(env);

  /* init */
  op = option_parser_new("[option ...] PNG_file [GFF3_file]",
                         "Create PNG representations of GFF3 annotation files.",
                         env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* -pipe */
  option = option_new_bool("pipe", "use pipe mode (i.e., show all gff3 "
                           "features on stdout)", &arguments->pipe, false, env);
  option_parser_add_option(op, option, env);

  /* -force */
  option = option_new_bool("force", "force writing to output file", &force,
                           false, env);
  option_parser_add_option(op, option, env);

  /* -seqid */
  option = option_new_string("seqid", "sequence region identifier"
                                      " (default: first one in file)",
                            arguments->seqid,
                            NULL, env);
  option_parser_add_option(op, option, env);
  option_hide_default(option);

  /* -start */
  option = option_new_ulong_min("start", "start position"
                                         " (default: region start)",
                            &arguments->start,
                            UNDEF_ULONG, 1, env);
  option_parser_add_option(op, option, env);
  option_hide_default(option);

  /* -end */
  option2 = option_new_ulong("end", "end position (default: region end)",
                            &arguments->end,
                            UNDEF_ULONG, env);
  option_parser_add_option(op, option2, env);
  /* -start and -end must be given together */
  option_imply(option, option2 ,env);
  option_imply(option2, option ,env);
  option_hide_default(option2);

  /* -width */
  option = option_new_uint_min("width", "target image width",
                            &arguments->width,
                            800, 1, env);
  option_parser_add_option(op, option, env);

  /* set contact mailaddress */
  option_parser_set_mailaddress(op, "<ssteinbiss@zbh.uni-hamburg.de>");

  /* parse options */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 2, env);

  if (oprval == OPTIONPARSER_OK && !force && file_exists(argv[*parsed_args])) {
    env_error_set(env, "file \"%s\" exists already. use option -force to "
                  "overwrite", argv[*parsed_args]);
    oprval = OPTIONPARSER_ERROR;
  }
  
  /* free */
  option_parser_delete(op, env);

  return oprval;
}

int gt_view(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream = NULL,
               *gff3_out_stream = NULL,
               *feature_stream = NULL;
  Gff3_view_arguments arguments;
  GenomeNode *gn = NULL;
  FeatureIndex *features = NULL;
  int parsed_args, had_err=0;
  char *seqid;
  Range qry_range, sequence_region_range;
  Array *results = NULL;
  Config *cfg;
  Str *config_file;
  char *prog;
  Splitter *splitter;

  env_error_check(env);

  /* option parsing */
  arguments.seqid = str_new(env);
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
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args + 1],
                                             arguments.verbose &&
                                             NULL, env);

  /* create gff3 output stream if -pipe was used */
  if (arguments.pipe)
    gff3_out_stream = gff3_out_stream_new(gff3_in_stream, NULL, env);

  /* create feature index and stream */
  features = feature_index_new(env);
  feature_stream = feature_stream_new((arguments.pipe ?
                                         gff3_out_stream :
                                         gff3_in_stream),
                                      features, env);

  /* check for correct order: range end < range start */
  if (!had_err &&
      arguments.start != UNDEF_ULONG &&
      arguments.end != UNDEF_ULONG &&
      !(arguments.start < arguments.end))
  {
    env_error_set(env, "start of query range (%lu) must be before "
                       "end of query range (%lu)",
                       arguments.start, arguments.end);
    had_err = -1;
  }

  if (!had_err)
  {
    /* pull the features through the stream and free them afterwards */
    while (!(had_err = genome_stream_next_tree(feature_stream, &gn, env)) && gn)
    {
      genome_node_rec_delete(gn, env);
    }
  }

  /* if seqid is empty, take first one added to index */
  if (!had_err && strcmp(str_get(arguments.seqid),"") == 0)
  {
    seqid = feature_index_get_first_seqid(features);
    if (seqid == NULL)
    {
      env_error_set(env, "GFF input file must contain a sequence region!");
      had_err = -1;
    }
  }
  else if (!had_err && !feature_index_has_seqid(features,
                                                str_get(arguments.seqid),
                                                env)) 
  {
    env_error_set(env, "sequence region '%s' does not exist in GFF input file",
                  str_get(arguments.seqid));
    had_err = -1;
  }
  else if (!had_err)
  {
    seqid = str_get(arguments.seqid);
  }

  results = array_new(sizeof (GenomeNode*), env);
  if (!had_err)
  {
    sequence_region_range = feature_index_get_range_for_seqid(features, seqid);
    qry_range.start = (arguments.start == UNDEF_ULONG ?
                         sequence_region_range.start :
                         arguments.start);
    qry_range.end   = (arguments.end == UNDEF_ULONG ?
                         sequence_region_range.end :
                         arguments.end);

    feature_index_get_features_for_range(features,
                                         results,
                                         seqid,
                                         qry_range,
                                         env);
  }

  if (!had_err)
  {
    if (arguments.verbose)
      fprintf(stderr, "# of results: %lu\n", array_size(results));

    /* find and load configuration file */
    prog = cstr_dup(argv[0], env); /* create modifiable copy for splitter */
    /* XXX: remove the ugly splitter stuff */
    splitter = splitter_new(env);
    splitter_split(splitter, prog, strlen(prog), ' ', env);
    config_file = gtdata_get_path(splitter_get_token(splitter, 0), env);
    splitter_delete(splitter, env);
    env_ma_free(prog, env);
    str_append_cstr(config_file, "/config/view.lua", env);
    cfg = config_new(env, arguments.verbose);
    assert(cfg);
    if (file_exists(str_get(config_file)))
      config_load_file(cfg, config_file, env);

    /* create and write image file */
    Diagram* d = diagram_new(results, qry_range, cfg, env);
    Render* r = render_new(cfg, env);
    render_to_png(r, d, (char*) argv[parsed_args], arguments.width, env);

    render_delete(r, env);
    config_delete(cfg, env);
    str_delete(config_file, env);
    diagram_delete(d, env);
  }

  /* free */
  str_delete(arguments.seqid,env);
  array_delete(results, env);
  feature_index_delete(features, env);
  genome_stream_delete(feature_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  genome_stream_delete(gff3_out_stream, env);

  return had_err;
}
