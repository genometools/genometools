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
#include "libgtext/cds_stream.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/seqid2file.h"

#define GT_CDS_SOURCE_TAG "gt cds"

typedef struct {
  bool verbose;
  Str *seqfile,
      *regionmapping;
} CDS_arguments;

static OPrval parse_options(int *parsed_args, CDS_arguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] GFF3_file", "Add CDS features to exon "
                         "features given in GFF3_file.", env);

  /* -seqfile and -regionmapping */
  seqid2file_options(op, arguments->seqfile, arguments->regionmapping, env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_cds(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream, *cds_stream = NULL, *gff3_out_stream = NULL;
  GenomeNode *gn;
  CDS_arguments arguments;
  RegionMapping *regionmapping;
  int parsed_args, had_err = 0;
  env_error_check(env);

  /* option parsing */
  arguments.seqfile = str_new();
  arguments.regionmapping = str_new();

  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(arguments.regionmapping);
      str_delete(arguments.seqfile);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(arguments.regionmapping);
      str_delete(arguments.seqfile);
      return 0;
  }

  /* create gff3 input stream */
  assert(parsed_args < argc);
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose, env);

  /* create region mapping */
  regionmapping = seqid2file_regionmapping_new(arguments.seqfile,
                                               arguments.regionmapping,
                                               env_error(env));
  if (!regionmapping)
    had_err = -1;

  /* create CDS stream */
  if (!had_err) {
    cds_stream = cds_stream_new(gff3_in_stream, regionmapping,
                                GT_CDS_SOURCE_TAG, env);
    if (!cds_stream)
      had_err = -1;
  }

  /* create gff3 output stream */
  /* XXX: replace NULL with proper outfile */
  if (!had_err)
    gff3_out_stream = gff3_out_stream_new(cds_stream, NULL, env);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn);
  }

  /* free */
  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(cds_stream);
  genome_stream_delete(gff3_in_stream);
  str_delete(arguments.regionmapping);
  str_delete(arguments.seqfile);

  return had_err;
}
