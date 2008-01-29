/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/outputfile.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtext/uniq_stream.h"
#include "tools/gt_uniq.h"

typedef struct {
  bool verbose;
  GenFile *outfp;
} UniqArguments;

static OPrval parse_options(int *parsed_args, UniqArguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  OutputFileInfo *ofi;
  Option *option;
  OPrval oprval;
  error_check(err);

  /* init */
  op = option_parser_new("[option ...] [GFF3_file]", "Filter out repeated "
                         "features in a sorted GFF3_file.");
  ofi = outputfileinfo_new();

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, ofi);

  /* parse options */
  option_parser_set_max_args(op, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);

  /* free */
  outputfileinfo_delete(ofi);
  option_parser_delete(op);

  return oprval;
}

int gt_uniq(int argc, const char **argv, Error *err)
{
  GenomeStream *gff3_in_stream,
               *uniq_stream = NULL,
               *gff3_out_stream = NULL;
  UniqArguments arguments;
  GenomeNode *gn;
  int parsed_args, had_err;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose &&
                                             arguments.outfp);

  /* create uniq stream */
  uniq_stream = uniq_stream_new(gff3_in_stream);

  /* create gff3 output stream */
  gff3_out_stream = gff3_out_stream_new(uniq_stream, arguments.outfp);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, err)) && gn)
    genome_node_rec_delete(gn);

  /* free */
  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(uniq_stream);
  genome_stream_delete(gff3_in_stream);
  genfile_close(arguments.outfp);

  return had_err;
}
