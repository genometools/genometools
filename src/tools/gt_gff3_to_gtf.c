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
#include "libgtcore/versionfunc.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gtf_out_stream.h"
#include "tools/gt_gff3_to_gtf.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Error *err)
{
  OptionParser *op;
  OPrval oprval;
  op = option_parser_new("[GFF3_file ...]", "Parse GFF3 file(s) and show it as "
                         "GTF2.2.");
  /* parse */
  option_parser_set_max_args(op, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_gff3_to_gtf(int argc, const char **argv, Error *err)
{
  GenomeStream *gff3_in_stream = NULL, *gtf_out_stream = NULL;
  GenomeNode *gn;
  int parsed_args, had_err = 0;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               false, false);

  if (!gff3_in_stream)
    had_err = -1;

  if (!had_err) {
    /* create a gtf output stream */
    gtf_out_stream = gtf_out_stream_new(gff3_in_stream, NULL);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = genome_stream_next_tree(gtf_out_stream, &gn, err)) &&
           gn) {
      genome_node_rec_delete(gn);
    }
  }

  /* free */
  genome_stream_delete(gff3_in_stream);
  genome_stream_delete(gtf_out_stream);

  return had_err;
}
