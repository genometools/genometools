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
#include "libgtext/regioncov_visitor.h"
#include "tools/gt_regioncov.h"

typedef struct {
  unsigned long max_feature_dist;
  bool verbose;
} RegionCovArguments;

static OPrval parse_options(int *parsed_args, RegionCovArguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] GFF3_file",
                         "Show which parts of the given sequence regions are "
                         "covered by features.");
  /* -maxfeaturedist */
  o = option_new_ulong("maxfeaturedist", "set the maximum distance two "
                       "features can have while still being in the same "
                       "``cluster''", &arguments->max_feature_dist, 0);
  option_parser_add_option(op, o);
  /* -v */
  o = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, o);
  /* parse */
  option_parser_set_min_max_args(op, 1, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_regioncov(int argc, const char **argv, Error *err)
{
  GenomeVisitor *regioncov_visitor;
  GenomeStream *gff3_in_stream;
  GenomeNode *gn;
  RegionCovArguments arguments;
  int parsed_args, had_err = 0;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      return 0;
  }

  /* create gff3 input stream */
  assert(parsed_args < argc);
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose);

  /* create region coverage visitor */
  regioncov_visitor = regioncov_visitor_new(arguments.max_feature_dist);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = genome_stream_next_tree(gff3_in_stream, &gn, err)) && gn) {
      had_err = genome_node_accept(gn, regioncov_visitor, err);
      genome_node_rec_delete(gn);
  }

  /* show region coverage */
  if (!had_err)
    regioncov_visitor_show_coverage(regioncov_visitor);

  /* free */
  genome_visitor_delete(regioncov_visitor);
  genome_stream_delete(gff3_in_stream);

  return had_err;
}
