/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

typedef struct {
  unsigned long max_feature_dist;
  bool verbose;
} RegionCovArguments;

static OPrval parse_options(int *parsed_args, RegionCovArguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] GFF3_file",
                         "Show which parts of the given sequence regions are "
                         "covered by features.", env);
  /* -maxfeaturedist */
  o = option_new_ulong("maxfeaturedist", "set the maximum distance two "
                       "features can have while still being in the same "
                       "``cluster''", &arguments->max_feature_dist, 0, env);
  option_parser_add_option(op, o, env);
  /* -v */
  o = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, o, env);
  /* parse */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_regioncov(int argc, const char **argv, Env *env)
{
  GenomeVisitor *regioncov_visitor;
  GenomeStream *gff3_in_stream;
  GenomeNode *gn;
  RegionCovArguments arguments;
  int parsed_args, had_err = 0;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      return 0;
  }

  /* create gff3 input stream */
  assert(parsed_args < argc);
  gff3_in_stream = gff3_in_stream_new_sorted(argv[parsed_args],
                                             arguments.verbose, env);

  /* create region coverage visitor */
  regioncov_visitor = regioncov_visitor_new(arguments.max_feature_dist, env);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = genome_stream_next_tree(gff3_in_stream, &gn, env)) && gn) {
      had_err = genome_node_accept(gn, regioncov_visitor, env);
      genome_node_rec_delete(gn, env);
  }

  /* show region coverage */
  if (!had_err)
    regioncov_visitor_show_coverage(regioncov_visitor, env);

  /* free */
  genome_visitor_delete(regioncov_visitor, env);
  genome_stream_delete(gff3_in_stream, env);

  return had_err;
}
