/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static int parse_options(FILE **outfp, int argc, char **argv)
{
  int parsed_args;
  OptionParser *op;
  Option *option;

  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Merge sorted GFF3 files in sorted fashion.");

  option = option_new_outputfile(outfp);
  option_parser_add_option(op, option);
  parsed_args = option_parser_parse(op, argc, argv, versionfunc);
  option_parser_free(op);
  return parsed_args;
}

int gt_merge(int argc, char *argv[])
{
  Genome_stream *gff3_in_stream,
                *merge_stream,
                *gff3_out_stream;
  Array *genome_streams = array_new(sizeof(Genome_stream*));
  Genome_node *gn;
  unsigned long i;
  int parsed_args;
  FILE *outfp;

  /* option parsing */
  parsed_args = parse_options(&outfp, argc, argv);

  /* XXX: check for multiple specification of '-' */

  /* create an gff3 input stream for each given file */
  if (parsed_args < argc) {
    /* we got files to open */
    for (i = parsed_args; i < argc; i++) {
      gff3_in_stream = gff3_in_stream_new_sorted(argv[i], 0);
      array_add(genome_streams, gff3_in_stream);
    }
   }
   else {
     /* use stdin */
     gff3_in_stream = gff3_in_stream_new_sorted(NULL, 0);
     array_add(genome_streams, gff3_in_stream);
   }

  /* create a merge stream */
  merge_stream = merge_stream_new(genome_streams);

  /* create a gff3 output stream */
  gff3_out_stream = gff3_out_stream_new(merge_stream, outfp);

  /* pull the features through the stream and free them afterwards */
  while ((gn = genome_stream_next_tree(gff3_out_stream, NULL)))
    genome_node_rec_free(gn);

  /* free */
  genome_stream_free(gff3_out_stream);
  genome_stream_free(merge_stream);
  for (i = 0; i < array_size(genome_streams); i++)
    genome_stream_free(*(Genome_stream**) array_get(genome_streams, i));
  array_free(genome_streams);
  if (outfp != stdout) xfclose(outfp);

  return EXIT_SUCCESS;
}
