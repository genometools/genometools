/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/extractfeat_stream.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/seqid2file.h"

typedef struct {
  bool join,
       translate,
       verbose;
  Str *type,
      *seqfile,
      *regionmapping;
} ExtractFeatArguments;

static OPrval parse_options(int *parsed_args, ExtractFeatArguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  error_check(err);
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

  /* -seqfile and -regionmapping */
  seqid2file_options(op, arguments->seqfile, arguments->regionmapping);

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, err);
  option_parser_delete(op);
  return oprval;
}

int gt_extractfeat(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream = NULL, *extractfeat_stream = NULL;
  GenomeNode *gn;
  GenomeFeatureType type;
  ExtractFeatArguments arguments;
  RegionMapping *regionmapping;
  int parsed_args, had_err = 0;
  env_error_check(env);

  /* option parsing */
  arguments.type = str_new();
  arguments.seqfile = str_new();
  arguments.regionmapping = str_new();
  switch (parse_options(&parsed_args, &arguments, argc, argv, env_error(env))) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(arguments.regionmapping);
      str_delete(arguments.seqfile);
      str_delete(arguments.type);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(arguments.regionmapping);
      str_delete(arguments.seqfile);
      str_delete(arguments.type);
      return 0;
  }

  /* determine type and make sure it is a valid one */
  if (genome_feature_type_get(&type, str_get(arguments.type)) == -1) {
    env_error_set(env, "\"%s\" is not a valid feature type",
                  str_get(arguments.type));
    had_err = -1;
  }

  if (!had_err) {
    /* create gff3 input stream */
    assert(parsed_args < argc);
    gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                               arguments.verbose);

    /* create region mapping */
    regionmapping = seqid2file_regionmapping_new(arguments.seqfile,
                                                 arguments.regionmapping,
                                                 env_error(env));
    if (!regionmapping)
      had_err = -1;
  }

  if (!had_err) {
    /* create extract feature stream */
    extractfeat_stream = extractfeat_stream_new(gff3_in_stream, regionmapping,
                                                type, arguments.join,
                                                arguments.translate);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = genome_stream_next_tree(extractfeat_stream, &gn,
                                               env_error(env))) && gn) {
      genome_node_rec_delete(gn);
    }
  }

  /* free */
  genome_stream_delete(extractfeat_stream);
  genome_stream_delete(gff3_in_stream);
  str_delete(arguments.regionmapping);
  str_delete(arguments.seqfile);
  str_delete(arguments.type);

  return had_err;
}
