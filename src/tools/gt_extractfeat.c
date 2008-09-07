/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/option.h"
#include "core/unused.h"
#include "extended/extract_feat_stream.h"
#include "extended/feature_type_factory_any.h"
#include "extended/gff3_in_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "tools/gt_extractfeat.h"

typedef struct {
  bool join,
       translate,
       verbose;
  GT_Str *typestr,
      *seqfile,
      *regionmapping;
  GT_GenomeFeatureType *type;
  GT_FeatureTypeFactory *feature_type_factory;
} ExtractFeatArguments;

static void* gt_extractfeat_arguments_new(void)
{
  ExtractFeatArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->typestr = gt_str_new();
  arguments->seqfile = gt_str_new();
  arguments->regionmapping = gt_str_new();
  arguments->feature_type_factory = gt_feature_type_factory_any_new();
  return arguments;
}

static void gt_extractfeat_arguments_delete(void *tool_arguments)
{
  ExtractFeatArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_feature_type_factory_delete(arguments->feature_type_factory);
  gt_str_delete(arguments->regionmapping);
  gt_str_delete(arguments->seqfile);
  gt_str_delete(arguments->typestr);
  ma_free(arguments);
}

static OptionParser* gt_extractfeat_option_parser_new(void *tool_arguments)
{
  ExtractFeatArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *option;
  assert(arguments);

  op = option_parser_new("[option ...] GFF3_file",
                         "Extract features given in GFF3_file from "
                         "sequence file.");

  /* -type */
  option = option_new_string("type", "set type of features to extract",
                             arguments->typestr, NULL);
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

  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  option_parser_set_min_max_args(op, 1, 1);

  return op;
}

static int gt_extractfeat_arguments_check(UNUSED int argc, void *tool_arguments,
                                          GT_Error *err)
{
  ExtractFeatArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  assert(arguments);

  /* determine type and make sure it is a valid one */
  if (!(arguments->type =
          gt_feature_type_factory_create_gft(arguments->feature_type_factory,
                                          gt_str_get(arguments->typestr)))) {
    gt_error_set(err, "\"%s\" is not a valid feature type",
              gt_str_get(arguments->typestr));
    had_err = -1;
  }

  return had_err;
}

static int gt_extractfeat_runner(UNUSED int argc, const char **argv,
                                 int parsed_args, void *tool_arguments,
                                 GT_Error *err)
{
  GenomeStream *gff3_in_stream = NULL, *extract_feat_stream = NULL;
  GT_GenomeNode *gn;
  ExtractFeatArguments *arguments = tool_arguments;
  RegionMapping *regionmapping;
  int had_err = 0;

  gt_error_check(err);
  assert(arguments);

  if (!had_err) {
    /* create gff3 input stream */
    gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                               arguments->verbose);
    gff3_in_stream_set_feature_type_factory(gff3_in_stream,
                                            arguments->feature_type_factory);

    /* create region mapping */
    regionmapping = seqid2file_regionmapping_new(arguments->seqfile,
                                                 arguments->regionmapping, err);
    if (!regionmapping)
      had_err = -1;
  }

  if (!had_err) {
    /* create extract feature stream */
    extract_feat_stream = extract_feat_stream_new(gff3_in_stream, regionmapping,
                                                  arguments->type,
                                                  arguments->join,
                                                  arguments->translate);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = genome_stream_next_tree(extract_feat_stream, &gn,
                                               err)) && gn) {
      gt_genome_node_rec_delete(gn);
    }
  }

  /* free */
  genome_stream_delete(extract_feat_stream);
  genome_stream_delete(gff3_in_stream);

  return had_err;
}

Tool* gt_extractfeat(void)
{
  return tool_new(gt_extractfeat_arguments_new,
                  gt_extractfeat_arguments_delete,
                  gt_extractfeat_option_parser_new,
                  gt_extractfeat_arguments_check,
                  gt_extractfeat_runner);
}
