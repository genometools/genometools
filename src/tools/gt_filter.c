/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool verbose;
  Str *seqid,
      *typefilter;
  unsigned long max_gene_length,
                max_gene_num;
  double min_gene_score;
  FILE *outfp;
} FilterArgumentss;

static OPrval parse_options(int *parsed_args, FilterArgumentss *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  env_error_check(env);

  op = option_parser_new("[option ...] [GFF3_file ...]", "Filter GFF3 files.",
                         env);

  /* -seqid */
  option = option_new_string("seqid", "seqid a feature must have to pass the "
                             "filter (excluding comments)", arguments->seqid,
                             NULL, env);
  option_parser_add_option(op, option, env);

  /* -typefilter */
  option = option_new_string("typefilter", "filter out all features of the "
                             "given type", arguments->typefilter, NULL, env);
  /* XXX */
  option_is_development_option(option);
  option_parser_add_option(op, option, env);

  /* -maxgenelength */
  option = option_new_ulong_min("maxgenelength", "the maximum length a gene "
                                "can have to pass the filter",
                                &arguments->max_gene_length, UNDEFULONG, 1,
                                env);
  option_parser_add_option(op, option, env);

  /* -maxgenenum */
  option = option_new_ulong("maxgenenum", "the maximum number of genes which "
                            "can pass the filter", &arguments->max_gene_num,
                            UNDEFULONG, env);
  option_parser_add_option(op, option, env);

  /* -mingenescore */
  option = option_new_double("mingenescore", "the minimum score a gene must "
                             "have to pass the filter",
                             &arguments->min_gene_score, UNDEFDOUBLE, env);
  option_parser_add_option(op, option, env);

  /* -o */
  option = option_new_outputfile(&arguments->outfp, env);
  option_parser_add_option(op, option, env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* parse */
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_filter(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream, *filter_stream, *gff3_out_stream;
  GenomeNode *gn;
  FilterArgumentss arguments;
  int parsed_args, has_err;

  /* option parsing */
  arguments.seqid = str_new(env);
  arguments.typefilter = str_new(env);
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(arguments.seqid, env);
      str_delete(arguments.typefilter, env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(arguments.seqid, env);
      str_delete(arguments.typefilter, env);
      return 0;
  }

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               arguments.verbose &&
                                               arguments.outfp != stdout, env);

  /* create a filter stream */
  filter_stream = filter_stream_new(gff3_in_stream, arguments.seqid,
                                    arguments.typefilter,
                                    arguments.max_gene_length,
                                    arguments.max_gene_num,
                                    arguments.min_gene_score, env);

  /* create a gff3 output stream */
  gff3_out_stream = gff3_out_stream_new(filter_stream, arguments.outfp, env);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_out_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn, env);
  }

  /* free */
  genome_stream_delete(gff3_out_stream, env);
  genome_stream_delete(filter_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  if (arguments.outfp != stdout)
    env_fa_xfclose(arguments.outfp, env);
  str_delete(arguments.seqid, env);
  str_delete(arguments.typefilter, env);

  return has_err;
}
