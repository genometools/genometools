/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(int *parsed_args, FILE **outfp, int argc,
                            char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Merge sorted GFF3 files in sorted fashion.", env);
  option = option_new_outputfile(outfp, env);
  option_parser_add_option(op, option, env);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_merge(int argc, char *argv[], Env *env)
{
  GenomeStream *gff3_in_stream,
                *merge_stream,
                *gff3_out_stream;
  Array *genome_streams;
  GenomeNode *gn;
  unsigned long i;
  int parsed_args, has_err;
  FILE *outfp;

  /* alloc */
  genome_streams = array_new(sizeof (GenomeStream*), env);

  /* option parsing */
  switch (parse_options(&parsed_args, &outfp, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* XXX: check for multiple specification of '-' */

  /* create an gff3 input stream for each given file */
  if (parsed_args < argc) {
    /* we got files to open */
    for (i = parsed_args; i < argc; i++) {
      gff3_in_stream = gff3_in_stream_new_sorted(argv[i], 0, env);
      array_add(genome_streams, gff3_in_stream, env);
    }
   }
   else {
     /* use stdin */
     gff3_in_stream = gff3_in_stream_new_sorted(NULL, 0, env);
     array_add(genome_streams, gff3_in_stream, env);
   }

  /* create a merge stream */
  merge_stream = merge_stream_new(genome_streams, env);

  /* create a gff3 output stream */
  gff3_out_stream = gff3_out_stream_new(merge_stream, outfp, env);

  /* pull the features through the stream and free them afterwards */
  while (!(has_err = genome_stream_next_tree(gff3_out_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn, env);
  }

  /* free */
  genome_stream_delete(gff3_out_stream, env);
  genome_stream_delete(merge_stream, env);
  for (i = 0; i < array_size(genome_streams); i++)
    genome_stream_delete(*(GenomeStream**) array_get(genome_streams, i), env);
  array_delete(genome_streams, env);
  if (outfp != stdout)
    xfclose(outfp);

  return has_err;
}
