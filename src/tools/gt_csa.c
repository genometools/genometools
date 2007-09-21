/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/csa_stream.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_out_stream.h"

#define DEFAULT_JOINLENGTH 300

typedef struct {
  bool verbose;
  unsigned long join_length;
  GenFile *outfp;
} CSAArguments;

static OPrval parse_options(int *parsed_args, CSAArguments *arguments,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  OutputFileInfo *ofi;
  Option *option;
  OPrval oprval;
  env_error_check(env);

  /* init */
  op = option_parser_new("[option ...] [GFF3_file]",
                         "Replace spliced alignments with computed consensus "
                         "spliced alignments.", env);
  ofi = outputfileinfo_new(env);

  /* -join-length */
  option = option_new_ulong("join-length", "set join length for the spliced "
                            "alignment clustering", &arguments->join_length,
                            DEFAULT_JOINLENGTH, env);
  option_parser_add_option(op, option, env);

  /* -v */
  option = option_new_verbose(&arguments->verbose, env);
  option_parser_add_option(op, option, env);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, ofi, env);

  /* parse options */
  oprval = option_parser_parse_max_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, env);

  /* free */
  outputfileinfo_delete(ofi, env);
  option_parser_delete(op, env);

  return oprval;
}

int gt_csa(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream,
               *csa_stream,
               *gff3_out_stream;
  GenomeNode *gn;
  CSAArguments arguments;
  int parsed_args, had_err;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* create the streams */
  gff3_in_stream  = gff3_in_stream_new_sorted(argv[parsed_args],
                                              arguments.verbose &&
                                              arguments.outfp, env);
  csa_stream      = csa_stream_new(gff3_in_stream, arguments.join_length, env);
  gff3_out_stream = gff3_out_stream_new(csa_stream, arguments.outfp, env);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = genome_stream_next_tree(gff3_out_stream, &gn, env)) &&
         gn) {
    genome_node_rec_delete(gn, env);
  }

  /* free */
  genome_stream_delete(gff3_out_stream, env);
  genome_stream_delete(csa_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  genfile_xclose(arguments.outfp, env);

  return had_err;
}
