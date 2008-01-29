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

#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtext/gtf_in_stream.h"
#include "tools/gt_gtf_to_gff3.h"

static OPrval parse_options(int *parsed_args, bool *be_tolerant, int argc,
                            const char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  op = option_parser_new("[gtf_file]",
                         "Parse GTF2.2 file and show it as GFF3.");
  /* -tolerant */
  option = option_new_bool("tolerant", "be tolerant when parsing the GTF file",
                           be_tolerant, false);
  option_parser_add_option(op, option);
  /* parse */
  option_parser_set_max_args(op, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_gtf_to_gff3(int argc, const char **argv, Error *err)
{
  GenomeStream *gtf_in_stream = NULL, *gff3_out_stream = NULL;
  GenomeNode *gn;
  int parsed_args, had_err = 0;
  bool be_tolerant;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &be_tolerant, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create a gtf input stream */
  gtf_in_stream = gtf_in_stream_new(argv[parsed_args], be_tolerant, err);
  if (!gtf_in_stream)
    had_err = -1;

  if (!had_err) {
    /* create a gff3 output stream */
    /* XXX: use proper genfile */
    gff3_out_stream = gff3_out_stream_new(gtf_in_stream, NULL);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, err)) &&
           gn) {
      genome_node_rec_delete(gn);
    }
  }

  /* free */
  genome_stream_delete(gff3_out_stream);
  genome_stream_delete(gtf_in_stream);

  return had_err;
}
